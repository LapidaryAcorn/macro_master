"""
"The KMV method" (DECISIONS Sec 31): KMV's estimation and discretization
procedure as read from their replication package (openICPSR 231376, local
only, never in git; paths and lines in DECISIONS Sec 31.1). Nothing here is
copied from that package: it is an independent Python implementation of the
procedure described there.

Paths below: E/ = EarningsProcess/EarningsEstimation/,
D/ = EarningsProcess/EarningsDiscretization/.

Contents
  KMV_SIM, KMV_WEIGHTS, KMV_BOUNDS, KMV_GUESS, KMV_TARGETS, KMV_PAPER_PARAMS
      the method's settings (process, design, moments, objective, box, start).
  dfls(...)
      derivative-free trust-region least squares (Gauss-Newton on linear
      interpolation models of the residual vector, in the spirit of DFO-GN,
      Cartis & Roberts 2017). It stands in for KMV's NEWUOA-H (Zhang, Conn &
      Scheinberg), which needs a Fortran compiler we do not have. Same class
      of method, same inputs: residual vector, rhobeg, rhoend, maxfun,
      sequential restarts.
  estimate_dfls(...)
      KMV's estimation: logistic transforms, MakeGuess start, ndfls restarts.
  kmv_grid / kmv_generator / kmv_discretize
      KMV's 3 x 11 chain: power-spaced grid whose width and spacing are fitted
      by DFLS to the 8 moments of the simulated chain.
  our_chain_kmv_design(...)
      our 5 x 15 discretization of the same (reset) process, scaled so that the
      chain's annual log variance matches the continuous process in KMV's
      ergodic design (the analogue of match_var_log).
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import expm
from scipy.special import expit, logit
from scipy.stats import norm

from .simulate import MOMENT_ORDER, PARAM_ORDER, model_moments, moments_from_panel

# ---------------------------------------------------------------- settings --

# process + design + moments of E/Simulate.f90, E/ComputeMoments.f90:
# reset (rho = 0), Euler decay, Bernoulli lambda*dt, dt = 0.25 quarter
# (4 steps per quarter = 16 per year), stationary start, 25-year burn-in
# (100 quarters), 5 annual observations, one cross-section, "d5" = year 5 -
# year 1 (a 4-year gap), nsim = 5000.
KMV_SIM = dict(n_workers=5_000, n_years_keep=5, steps_per_quarter=4,
               lifecycle=False, burn_in_years=25, reset_jumps=True,
               kmv_euler=True, cross_section=True, d5_gap=4)

# E/dfovec.f90 l. 126, 145, 162: weight 0.5 on both kurtoses
KMV_WEIGHTS = {m: (0.5 if m.startswith("kurt") else 1.0) for m in MOMENT_ORDER}

# E/SetParameters.f90 l. 65-69: lambda, sigma < 2, beta < 1 (logistic, so the
# lower limit is 0). For log-scale DE/TikTak the lower ends are ours.
KMV_UPPER = {"lambda1": 2.0, "lambda2": 2.0, "sigma1": 2.0, "sigma2": 2.0,
             "beta1": 1.0, "beta2": 1.0}
KMV_BOUNDS = {"lambda1": (0.01, 2.0), "beta1": (0.05, 1.0), "sigma1": (0.2, 2.0),
              "lambda2": (0.001, 2.0), "beta2": (1e-6, 1.0), "sigma2": (0.2, 2.0)}

# E/MakeGuess.f90 l. 9-19
KMV_GUESS = {"lambda1": 0.08, "beta1": 0.7, "sigma1": 1.6,
             "lambda2": 0.007, "beta2": 0.01, "sigma2": 1.6}

# E/SetParameters.f90 l. 49-59 (Table 2 prints 11.6 for kurt d5)
KMV_TARGETS = {"var_log_earns": 0.70, "var_d1": 0.23, "var_d5": 0.46,
               "kurt_d1": 17.8, "kurt_d5": 11.55, "frac_d1_lt_10": 0.54,
               "frac_d1_lt_20": 0.71, "frac_d1_lt_50": 0.86}

# E/earnings_estimation_paper/parameters.txt (beta = "delta" in their code)
KMV_PAPER_PARAMS = {"lambda1": 0.07987252454850691, "beta1": 0.761606393760105,
                    "sigma1": 1.73522438989472, "lambda2": 0.006560658609206301,
                    "beta2": 0.009391785115626633, "sigma2": 1.52585125038250}
# E/earnings_estimation_output/parameters.txt: KMV's own re-run
KMV_RERUN_PARAMS = {"lambda1": 0.07988762707211912, "beta1": 0.754396874946948,
                    "sigma1": 1.73405631833937, "lambda2": 0.006528117159265131,
                    "beta2": 0.009205963050942997, "sigma2": 1.52730372411936}
# Table 3 (rounded, quarterly)
KMV_TABLE3 = {"lambda1": 0.080, "beta1": 0.761, "sigma1": 1.74,
              "lambda2": 0.007, "beta2": 0.009, "sigma2": 1.53}

# KMV's estimation order of the parameter vector (E/Estimation.f90 l. 32-79)
KMV_PARAM_ORDER = ["lambda1", "lambda2", "sigma1", "sigma2", "beta1", "beta2"]

DFLS_EST = dict(ndfls=3, rhobeg=5.0, rhoend=1e-4)       # E/SetParameters l. 61-63
DFLS_DISC = dict(ndfls=4, rhobeg=5.0, rhoend=1e-4)      # D/SetParameters l. 59-61


def kmv_objective_value(m: dict, targets: dict, weights: dict | None = None) -> float:
    w = KMV_WEIGHTS if weights is None else weights
    return float(sum(w[k] * (m[k] / targets[k] - 1.0) ** 2 for k in MOMENT_ORDER))


def relabel(p: dict) -> dict:
    """Components are exchangeable under KMV's box (both lambdas < 2): call
    component 1 the one with the larger arrival rate."""
    if p["lambda2"] > p["lambda1"]:
        return {"lambda1": p["lambda2"], "beta1": p["beta2"], "sigma1": p["sigma2"],
                "lambda2": p["lambda1"], "beta2": p["beta1"], "sigma2": p["sigma1"]}
    return dict(p)


# -------------------------------------------------------------------- DFLS --

def _tr_gn_step(J: np.ndarray, r: np.ndarray, delta: float) -> np.ndarray:
    """min ||r + J s||^2 s.t. ||s|| <= delta (Levenberg-Marquardt / More-Sorensen
    on the small system)."""
    g = J.T @ r
    H = J.T @ J
    n = len(g)
    try:
        s = -np.linalg.lstsq(J, r, rcond=None)[0]
    except np.linalg.LinAlgError:
        s = -g
    if np.linalg.norm(s) <= delta:
        return s
    lo, hi = 0.0, max(1e-12, np.linalg.norm(g) / delta + np.linalg.norm(H, 2))
    for _ in range(200):                       # bisection on mu
        mu = 0.5 * (lo + hi)
        s = -np.linalg.solve(H + mu * np.eye(n), g)
        if np.linalg.norm(s) > delta:
            lo = mu
        else:
            hi = mu
        if hi - lo < 1e-12 * max(1.0, hi):
            break
    s = -np.linalg.solve(H + hi * np.eye(n), g)
    return s


def dfls(resid, x0: np.ndarray, rhobeg: float = 5.0, rhoend: float = 1e-4,
         maxfun: int | None = None, log: list | None = None,
         trace: list | None = None):
    """Minimize sum(resid(x)**2) without derivatives (one DFLS attempt).

    Linear interpolation model of the residual vector on n + 1 points, a
    Gauss-Newton trust-region step, Powell-style radius management (rho
    falls from rhobeg to rhoend, the trust radius delta >= rho). Returns
    (x_best, f_best, n_evals).
    """
    x0 = np.asarray(x0, float)
    n = len(x0)
    maxfun = maxfun or 500 * (n + 1)
    nf = 0

    def ev(x):
        nonlocal nf
        nf += 1
        r = np.asarray(resid(x), float)
        f = float(r @ r)
        if log is not None:
            log.append((x.copy(), f))
        return r, f

    Y = [x0.copy()] + [x0 + rhobeg * e for e in np.eye(n)]
    RF = [ev(y) for y in Y]
    R = [rf[0] for rf in RF]
    F = [rf[1] for rf in RF]
    rho = delta = rhobeg

    def reduce_rho(rho):
        if rho <= rhoend:
            return None
        if rho > 250 * rhoend:
            return rho / 10.0
        if rho > 16 * rhoend:
            return np.sqrt(rho * rhoend)
        return rhoend

    while nf < maxfun:
        k = int(np.argmin(F))
        xk, rk, fk = Y[k], R[k], F[k]
        others = [i for i in range(n + 1) if i != k]
        D = np.array([Y[i] - xk for i in others])            # n x n displacements
        dR = np.array([R[i] - rk for i in others])            # n x m
        try:
            J = np.linalg.solve(D, dR).T                       # m x n
            cond_ok = np.linalg.cond(D) < 1e8
        except np.linalg.LinAlgError:
            J, cond_ok = None, False
        dist = np.array([np.linalg.norm(Y[i] - xk) for i in others])

        if not cond_ok:                                       # geometry repair
            t = others[int(np.argmax(dist))]
            Dinv = np.linalg.pinv(D)
            gl = Dinv[:, others.index(t)]                      # Lagrange gradient of t
            d = gl / max(np.linalg.norm(gl), 1e-300)
            Y[t] = xk + delta * d
            R[t], F[t] = ev(Y[t])
            continue

        s = _tr_gn_step(J, rk, delta)
        ns = np.linalg.norm(s)
        if ns < 0.5 * rho:                                    # step too short
            delta = max(0.5 * delta, rho)
            if dist.max() > 2.0 * delta:                      # fix the model first
                t = others[int(np.argmax(dist))]
                gl = np.linalg.inv(D)[:, others.index(t)]
                Y[t] = xk + delta * gl / max(np.linalg.norm(gl), 1e-300)
                R[t], F[t] = ev(Y[t])
                continue
            new = reduce_rho(rho)
            if new is None:
                break
            rho = new
            delta = max(0.5 * delta, rho)
            continue

        xn = xk + s
        rn, fn = ev(xn)
        pred = fk - float((rk + J @ s) @ (rk + J @ s))
        ratio = (fk - fn) / pred if pred > 0 else -1.0
        # Powell's NEWUOA radius rule (newuob: ratio <= 0.1 / <= 0.7 / > 0.7)
        if ratio <= 0.1:
            delta = 0.5 * ns
        elif ratio <= 0.7:
            delta = max(0.5 * delta, ns)
        else:
            delta = max(0.5 * delta, 2.0 * ns)
        if delta <= 1.5 * rho:
            delta = rho

        # point to drop: Powell's criterion |l_t(xn)| * max(1, (dist_t/delta)^2)
        Dinv = np.linalg.inv(D)
        lag = (xn - xk) @ Dinv                                 # l_t(xn), t in others
        lag_k = 1.0 - lag.sum()
        if fn < fk:
            cand = list(range(n + 1))
            vals = []
            for i in cand:
                l_i = lag_k if i == k else lag[others.index(i)]
                di = np.linalg.norm(Y[i] - xn)
                vals.append(abs(l_i) * max(1.0, (di / delta) ** 2))
            t = cand[int(np.argmax(vals))]
        else:
            vals = [abs(lag[j]) * max(1.0, (dist[j] / delta) ** 2) for j in range(n)]
            t = others[int(np.argmax(vals))]
        Y[t], R[t], F[t] = xn, rn, fn

        if ratio <= 0.1 and delta <= rho:                     # failure at the floor
            kb = int(np.argmin(F))
            dd = [np.linalg.norm(y - Y[kb]) for y in Y]
            if max(dd) > 2.0 * rho:                           # model step: move the
                t = int(np.argmax(dd))                        # farthest point close
                ob = [i for i in range(n + 1) if i != kb]
                Db = np.array([Y[i] - Y[kb] for i in ob])
                gl = np.linalg.pinv(Db)[:, ob.index(t)]
                Y[t] = Y[kb] + rho * gl / max(np.linalg.norm(gl), 1e-300)
                R[t], F[t] = ev(Y[t])
                continue
            new = reduce_rho(rho)
            if new is None:
                break
            rho = new
            delta = max(0.5 * delta, rho)
        if trace is not None:
            trace.append((nf, rho, delta, float(min(F))))

    k = int(np.argmin(F))
    return Y[k].copy(), F[k], nf


# ----------------------------------------------------- KMV's estimation ------

def _to_x(p: dict, order=KMV_PARAM_ORDER) -> np.ndarray:
    return np.array([logit(p[k] / KMV_UPPER[k]) for k in order])


def _to_p(x: np.ndarray, order=KMV_PARAM_ORDER) -> dict:
    return {k: float(KMV_UPPER[k] * expit(v)) for k, v in zip(order, x)}


class KMVResiduals:
    """sqrt(w_m) * (model_m / target_m - 1), E/dfovec.f90 l. 126-185 (picklable)."""

    def __init__(self, targets: dict, sim_kwargs: dict | None = None,
                 weights: dict | None = None):
        self.t = dict(targets)
        self.sim = dict(KMV_SIM if sim_kwargs is None else sim_kwargs)
        self.w = KMV_WEIGHTS if weights is None else weights

    def moments(self, p: dict) -> dict:
        return model_moments(p, **self.sim)

    def __call__(self, x: np.ndarray) -> np.ndarray:
        m = self.moments(_to_p(x))
        return np.array([np.sqrt(self.w[k]) * (m[k] / self.t[k] - 1.0) for k in MOMENT_ORDER])


def estimate_dfls(targets: dict, sim_kwargs: dict | None = None,
                  weights: dict | None = None, guess: dict | None = None,
                  ndfls: int = 3, rhobeg: float = 5.0, rhoend: float = 1e-4,
                  maxfun: int | None = None, keep_log: bool = False):
    """KMV's estimation (E/Estimation.f90 l. 98-114): ndfls sequential DFLS
    attempts, each restarting from the previous end point, in logistic
    coordinates. Returns (params, info)."""
    res = KMVResiduals(targets, sim_kwargs, weights)
    x = _to_x(guess or KMV_GUESS)
    log = [] if keep_log else None
    attempts = []
    for j in range(ndfls):
        x, f, nf = dfls(res, x, rhobeg=rhobeg, rhoend=rhoend, maxfun=maxfun, log=log)
        attempts.append({"attempt": j + 1, "f": f, "n_evals": nf, "params": _to_p(x)})
    p = _to_p(x)
    info = {"f": attempts[-1]["f"], "attempts": attempts,
            "n_evals": int(sum(a["n_evals"] for a in attempts))}
    if keep_log:
        info["log"] = log
    return p, info


# ------------------------------------------------- KMV's discretization ------

def kmv_grid(n: int, width: float, gridpar: float) -> np.ndarray:
    """D/TransitionMatrix.f90 l. 36-53 + PowerSpacedGrid (D/Procedures.f90
    l. 448-483): upper half 0 .. width/2 spaced as x^(1/gridpar), mirrored."""
    m = (n + 1) // 2
    x = np.linspace(0.0, 1.0, m)
    upper = (width / 2.0) * x ** (1.0 / gridpar)
    g = np.empty(n)
    g[m - 1:] = upper
    g[:m - 1] = -upper[1:][::-1]
    return g


def _find_lin_prob1(xi: float, x: np.ndarray) -> tuple[int, int]:
    """Indices (0-based) of the bracket of xi in x, as FindLinProb1
    (D/Procedures.f90 l. 556-597)."""
    n = len(x)
    below = np.where(xi > x)[0]
    if xi <= x[0]:
        return 0, 1
    locl = int(below[np.argmax(x[below])])        # MAXLOC(x, MASK = xi > x)
    if locl >= n - 1:
        return n - 2, n - 1
    return locl, locl + 1


def kmv_generator(grid: np.ndarray, lam: float, beta: float, sigma: float,
                  rho: float = 0.0) -> np.ndarray:
    """Continuous-time generator of one component, D/TransitionMatrix.f90:
    jumps l. 75-98/100-123 (mid-point bins, N(rho*origin, sigma^2), rows
    renormalized), drift l. 125-173 (DriftPointsVersion 2, deltaforapprox 1)."""
    n = len(grid)
    dg = np.diff(grid)
    J = np.empty((n, n))
    for i in range(n):
        mu = rho * grid[i]
        for j in range(n):
            if j == 0:
                J[i, j] = norm.cdf(grid[j] + 0.5 * dg[j], mu, sigma)
            elif j < n - 1:
                J[i, j] = (norm.cdf(grid[j] + 0.5 * dg[j], mu, sigma)
                           - norm.cdf(grid[j] - 0.5 * dg[j - 1], mu, sigma))
            else:
                J[i, j] = 1.0 - norm.cdf(grid[j] - 0.5 * dg[j - 1], mu, sigma)
        J[i] /= J[i].sum()
    Q = lam * (J - np.eye(n))
    Dm = np.zeros((n, n))
    for i in range(n):
        y = grid[i]
        if y == 0.0:
            continue
        lo, hi = _find_lin_prob1((1.0 - beta) * y, grid)
        Dm[i, i] = -1.0
        if y < 0.0:
            Dm[i, i] += (grid[hi] - (1.0 - beta) * y) / (grid[hi] - y)
            Dm[i, hi] = (-y + (1.0 - beta) * y) / (grid[hi] - y)
        else:
            Dm[i, lo] = (y - (1.0 - beta) * y) / (y - grid[lo])
            Dm[i, i] += (-grid[lo] + (1.0 - beta) * y) / (y - grid[lo])
    return Q + Dm


def generator_stationary(Q: np.ndarray) -> np.ndarray:
    n = Q.shape[0]
    A = np.vstack([Q.T, np.ones(n)])
    b = np.zeros(n + 1)
    b[-1] = 1.0
    pi = np.linalg.lstsq(A, b, rcond=None)[0]
    pi = np.clip(pi, 0.0, None)
    return pi / pi.sum()


def kmv_chain_panel(g1, Q1, pi1, g2, Q2, pi2, nsim: int = 50_000, seed: int = 7755,
                    dt: float = 0.25) -> np.ndarray:
    """Log annual earnings (nsim x 5) of the chain as D/Simulate.f90 simulates it:
    Euler transition I + dt*Q per step (l. 21-22), first step drawn from the
    ergodic distribution (l. 56-57), 81 steps, year k = steps 16(k-1)+1 .. 16k
    (l. 68-71)."""
    rng = np.random.default_rng(seed)
    T = int(np.floor(20.0 / dt)) + 1
    u1 = rng.random((nsim, T))
    u2 = rng.random((nsim, T))
    out = []
    for g, Q, pi, u in ((g1, Q1, pi1, u1), (g2, Q2, pi2, u2)):
        Pdt = np.eye(len(g)) + dt * Q
        if (Pdt < -1e-12).any():
            raise ValueError("I + dt*Q has negative entries (KMV's code would pause here)")
        cum = np.cumsum(np.clip(Pdt, 0.0, None), axis=1)
        cum /= cum[:, -1:]
        cpi = np.cumsum(pi)
        cpi /= cpi[-1]
        s = np.empty((nsim, T), dtype=np.int64)
        s[:, 0] = np.searchsorted(cpi, u[:, 0], side="left")
        for t in range(1, T):
            s[:, t] = (u[:, t:t + 1] > cum[s[:, t - 1]]).sum(axis=1)
        out.append(g[np.minimum(s, len(g) - 1)])
    lev = np.exp(out[0] + out[1])
    spy = int(round(4.0 / dt))
    ann = np.add.reduceat(lev[:, :5 * spy], np.arange(0, 5 * spy, spy), axis=1)
    return np.log(ann)


def kmv_chain_from_grid_params(params: dict, widths, gridpars, n1: int = 3,
                               n2: int = 11) -> dict:
    out = {}
    for c, n, w, k in (("1", n1, widths[0], gridpars[0]), ("2", n2, widths[1], gridpars[1])):
        g = kmv_grid(n, w, k)
        Q = kmv_generator(g, params[f"lambda{c}"], params[f"beta{c}"], params[f"sigma{c}"])
        pi = generator_stationary(Q)
        out[c] = {"grid": g, "Q": Q, "P": expm(Q), "pi": pi}
    return out


class _DiscResiduals:
    def __init__(self, params, targets, n1, n2, nsim, seed, weights):
        self.p, self.t, self.n1, self.n2 = params, targets, n1, n2
        self.nsim, self.seed = nsim, seed
        self.w = KMV_WEIGHTS if weights is None else weights
        s1, s2 = params["sigma1"], params["sigma2"]
        self.lo = np.array([s1, s2, 0.3, 0.3])               # D/MakeGuess.f90 l. 8-20
        self.hi = np.array([10 * s1, 10 * s2, 1.0, 1.0])

    def unpack(self, x):
        v = self.lo + (self.hi - self.lo) * expit(x)
        return v[:2], v[2:]

    def moments(self, x):
        widths, gps = self.unpack(x)
        c = kmv_chain_from_grid_params(self.p, widths, gps, self.n1, self.n2)
        y = kmv_chain_panel(c["1"]["grid"], c["1"]["Q"], c["1"]["pi"],
                            c["2"]["grid"], c["2"]["Q"], c["2"]["pi"],
                            nsim=self.nsim, seed=self.seed)
        return moments_from_panel(y, cross_section=True, d5_gap=4)

    def __call__(self, x):
        try:
            m = self.moments(x)
        except ValueError:
            return np.full(len(MOMENT_ORDER), 10.0)
        return np.array([np.sqrt(self.w[k]) * (m[k] / self.t[k] - 1.0) for k in MOMENT_ORDER])


def kmv_discretize(params: dict, targets: dict, n1: int = 3, n2: int = 11,
                   nsim: int = 50_000, seed: int = 7755, ndfls: int = 4,
                   rhobeg: float = 5.0, rhoend: float = 1e-4, weights: dict | None = None):
    """KMV's discretization (D/Estimation.f90): fit (width1, width2, gridpar1,
    gridpar2) by DFLS so the simulated chain's 8 moments match `targets` (the
    code uses the data targets, D/SetParameters.f90 l. 47-57). Returns
    (disc, info); `disc` has the layout of discretize.discretize_process, so
    discretize.export_chain / diagnose work on it."""
    from .discretize import _combine
    r = _DiscResiduals(params, targets, n1, n2, nsim, seed, weights)
    s1, s2 = params["sigma1"], params["sigma2"]
    g0 = np.array([2 * s1, 2 * s2, 0.7, 0.7])
    x = logit((g0 - r.lo) / (r.hi - r.lo))
    attempts = []
    for j in range(ndfls):
        x, f, nf = dfls(r, x, rhobeg=rhobeg, rhoend=rhoend)
        attempts.append({"attempt": j + 1, "f": f, "n_evals": nf})
    widths, gps = r.unpack(x)
    disc = kmv_chain_from_grid_params(params, widths, gps, n1, n2)
    disc = _combine(disc)
    disc["reset_jumps"] = True
    info = {"widths": [float(v) for v in widths], "gridpars": [float(v) for v in gps],
            "f": attempts[-1]["f"], "attempts": attempts,
            "moments": r.moments(x)}
    return disc, info


def chain_ergodic_var(disc: dict) -> dict:
    c = disc["combined"]
    z, pi = c["zgrid"], c["pi"]
    out = {"total": float(pi @ (z - pi @ z) ** 2)}
    for k in ("1", "2"):
        g, p = disc[k]["grid"], disc[k]["pi"]
        out[f"z{k}"] = float(p @ (g - p @ g) ** 2)
    return out


def our_chain_kmv_design(params: dict, n1: int = 5, n2: int = 15,
                         sim_kwargs: dict | None = None, nsim_chain: int = 50_000,
                         reset_jumps: bool = True):
    """Our 5 x 15 scheme (discretize.discretize_process: conditional-mean grid,
    analytic variance rescale) on the estimated process, then one global scale
    so that the chain's year-1 cross-sectional var(log annual earnings),
    simulated from its ergodic distribution, equals the continuous process's
    in `sim_kwargs` (default KMV_SIM with 200k workers). The analogue of
    match_var_log for an ergodic design."""
    from .discretize import discretize_process, simulate_discrete_panel, _combine
    sk = dict(KMV_SIM, n_workers=200_000)
    sk.update(sim_kwargs or {})
    disc = discretize_process(params, n1=n1, n2=n2, reset_jumps=reset_jumps)
    target = model_moments(params, **sk)["var_log_earns"]
    y = simulate_discrete_panel(disc, n_workers=nsim_chain, n_years_keep=5, lifecycle=False)
    got = float(np.var(y[:, 0]))
    s = np.sqrt(target / got)
    for k in ("1", "2"):
        disc[k]["grid"] = disc[k]["grid"] * s
    disc = _combine({k: disc[k] for k in ("1", "2")})
    if reset_jumps:
        disc["reset_jumps"] = True
    return disc, {"scale": float(s), "target_var_log": target, "chain_var_log_before": got}
