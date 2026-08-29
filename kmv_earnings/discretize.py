"""
Discretization of the estimated jump-drift process into a finite-state Markov
chain (the "Model Discretized" column of KMV Table 2), plus export of the chain
for downstream models (KMV HANK codes / Auclert et al. sequence-space toolkit).

Scheme, per component j:
  - bin edges power-spaced on [-zmax_mult*sigma_j, +zmax_mult*sigma_j]
    (outermost edges extended to +-inf); grid POINTS are the conditional means
    of N(0, sigma_j^2) within each bin, so the size of a discrete move matches
    the average jump landing in that cell;
  - drift -beta_j*z: upwind finite differences;
  - ADDITIVE jumps at rate lambda_j: landing probabilities integrate
    N(current state, sigma_j^2) over the bins;
  - generator Q_j -> quarterly transition matrix expm(Q_j);
  - VARIANCE CORRECTION (see below), then the combined chain is the Kronecker
    product over the two independent components.

Why the variance correction
---------------------------
The upwind scheme is numerically diffusive: probability mass that should drift
smoothly toward zero instead hops a whole grid point at a time, which inflates
the chain's dispersion. This bias is systematic - it does NOT vanish as the
grid is refined. At KMV's parameters the raw chain has ergodic variance ~1.27
against a theoretical lambda*sigma^2/(2*beta) total of ~1.07 (~19% too high),
essentially unchanged from a 33-state to a 465-state grid.

The fix (`rescale_variance=True`, default) rescales each component's grid by
sqrt(v_theory / v_chain), where v_theory = lambda_j*sigma_j^2/(2*beta_j) is the
exact stationary variance of the jump-drift component and v_chain is the
chain's. This is exact, analytic, and cheap. Crucially, KURTOSIS IS
SCALE-INVARIANT, so the correction fixes the variances without disturbing the
tail behaviour recovered by the grid.

Grid size and the kurtosis
--------------------------
Kurtosis needs grid points far enough out to represent rare large jumps. KMV's
3 x 11 = 33 states are too coarse for that here; the default is 5 x 15 = 75.
Measured against the continuous simulation at KMV's parameters
(0.69 / 0.23 / 0.50 / 14.60 / 10.76 / 0.55 / 0.66 / 0.84):

    (n1, n2)  states   var_log  var_d1  var_d5  kurt_d1  kurt_d5
    (3, 11)     33      0.63     0.19    0.45    10.91     9.45
    (5, 15)     75      0.64     0.20    0.45    13.87    10.24
    (9, 21)    189      0.64     0.20    0.46    14.74    10.64
    (11, 25)   275      0.64     0.20    0.46    14.78    10.69

Raising the grid past ~(9, 21) buys nothing. Larger chains are more faithful
but slow the step-2 solve, so the default trades off at 75 states; pass
n1/n2 explicitly if step 2 can afford more (or must have fewer).

Known residual gaps (honest limitations)
----------------------------------------
1. Variances land ~7% below the continuous process (0.64 vs 0.69). The
   analytic correction targets the ERGODIC variance, while the moments are
   computed on a 36-year lifecycle panel, and the chain's transient mixing
   differs slightly from the continuous process. Set `match_var_log=True` to
   calibrate a single extra scale factor by simulation so that the chain
   reproduces the continuous var(log annual earnings) - this costs two
   simulations and closes the remaining gap.
2. "Frac 1yr change < 10%" comes out too HIGH (~0.62 vs 0.55). This is
   intrinsic to discretization: in many quarters the chain does not change
   state at all, producing annual changes of exactly zero and piling up mass
   at tiny changes. No grid refinement removes it.
3. For an exact reproduction of KMV's published discretized column, substitute
   the routine from their replication package (appendix D.1); nothing else in
   this pipeline depends on how P is built.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import expm
from scipy.stats import norm

from .simulate import moments_from_panel, model_moments


def component_grid_and_edges(n: int, sigma: float, zmax_mult: float = 3.5,
                             curv: float = 1.5):
    """Bin edges power-spaced in [-zmax_mult*sigma, zmax_mult*sigma] (outer
    edges -> +-inf); grid points = conditional means of N(0, sigma^2) per bin."""
    x = np.linspace(-1.0, 1.0, n + 1)
    edges = np.sign(x) * np.abs(x) ** curv * zmax_mult * sigma
    edges[0], edges[-1] = -np.inf, np.inf
    pts = np.empty(n)
    for k in range(n):
        a, b = edges[k], edges[k + 1]
        mass = norm.cdf(b, scale=sigma) - norm.cdf(a, scale=sigma)
        pa = norm.pdf(a / sigma) if np.isfinite(a) else 0.0
        pb = norm.pdf(b / sigma) if np.isfinite(b) else 0.0
        pts[k] = sigma * (pa - pb) / max(mass, 1e-12)
    return pts, edges


def component_generator(grid: np.ndarray, edges: np.ndarray, lam: float,
                        beta: float, sigma: float) -> np.ndarray:
    n = len(grid)
    Q = np.zeros((n, n))
    for k in range(n):                       # drift toward zero (upwind)
        g = grid[k]
        if g > 0 and k > 0:
            r = beta * g / (grid[k] - grid[k - 1])
            Q[k, k - 1] += r
            Q[k, k] -= r
        elif g < 0 and k < n - 1:
            r = beta * (-g) / (grid[k + 1] - grid[k])
            Q[k, k + 1] += r
            Q[k, k] -= r
    for k in range(n):                       # additive jumps
        w = np.diff(norm.cdf(edges, loc=grid[k], scale=sigma))
        Q[k, :] += lam * w
        Q[k, k] -= lam
    return Q


def stationary_dist(P: np.ndarray) -> np.ndarray:
    vals, vecs = np.linalg.eig(P.T)
    k = int(np.argmin(np.abs(vals - 1.0)))
    pi = np.abs(np.real(vecs[:, k]))
    return pi / pi.sum()


def theoretical_stationary_var(lam: float, beta: float, sigma: float) -> float:
    """Exact stationary variance of dz = -beta z dt + eps dN: lam*sigma^2/(2 beta)."""
    return lam * sigma ** 2 / (2.0 * beta)


def _combine(out: dict) -> dict:
    g1, g2 = out["1"]["grid"], out["2"]["grid"]
    zgrid = (g1[:, None] + g2[None, :]).ravel()
    P = np.kron(out["1"]["P"], out["2"]["P"])
    pi = np.kron(out["1"]["pi"], out["2"]["pi"])
    e = np.exp(zgrid)
    e = e / (pi @ e)                          # normalise mean earnings to 1
    out["combined"] = {"zgrid": zgrid, "egrid": e, "P": P, "pi": pi}
    return out


def discretize_process(params: dict, n1: int = 5, n2: int = 15,
                       zmax_mult: float = 3.5, curv: float = 1.5,
                       rescale_variance: bool = True,
                       match_var_log: bool = False,
                       sim_kwargs: dict | None = None):
    """
    Discretize both components and build the combined chain.

    n1 / n2      grid points for the transitory / persistent component
                 (default 5 x 15 = 75 states; KMV used 3 x 11 = 33, which is
                 too coarse to carry the kurtosis - see module docstring).
    rescale_variance  correct the upwind scheme's variance inflation by
                 matching each component's exact stationary variance (default
                 True; kurtosis is unaffected because it is scale-invariant).
    match_var_log     additionally calibrate ONE global scale factor by
                 simulation so the chain reproduces the continuous process's
                 var(log annual earnings). Costs two simulations.
    """
    specs = {
        "1": (params["lambda1"], params["beta1"], params["sigma1"], n1),
        "2": (params["lambda2"], params["beta2"], params["sigma2"], n2),
    }
    out = {}
    for name, (lam, beta, sig, n) in specs.items():
        grid, edges = component_grid_and_edges(n, sig, zmax_mult, curv)
        Q = component_generator(grid, edges, lam, beta, sig)
        P = expm(Q)                           # quarterly transition matrix
        pi = stationary_dist(P)
        if rescale_variance:
            v_chain = float(pi @ (grid - pi @ grid) ** 2)
            v_theory = theoretical_stationary_var(lam, beta, sig)
            if v_chain > 0:
                grid = grid * np.sqrt(v_theory / v_chain)
        out[name] = {"grid": grid, "edges": edges, "Q": Q, "P": P, "pi": pi}

    if match_var_log:
        sk = dict(n_workers=30_000, n_years_keep=36)
        sk.update(sim_kwargs or {})
        target = model_moments(params, **sk)["var_log_earns"]
        got = moments_from_panel(
            simulate_discrete_panel(_combine(dict(out)), n_workers=sk["n_workers"],
                                    n_years_keep=sk["n_years_keep"])
        )["var_log_earns"]
        if got > 0:
            s = np.sqrt(target / got)
            for name in ("1", "2"):
                out[name]["grid"] = out[name]["grid"] * s

    return _combine(out)


def simulate_discrete_panel(disc: dict, n_workers: int = 50_000,
                            n_years_keep: int = 36, lifecycle: bool = True,
                            seed: int = 4321) -> np.ndarray:
    """
    Simulate the chain at its native quarterly frequency (the frequency the
    HANK model uses), aggregate the earnings flow exp(z) to annual, return log
    annual earnings. lifecycle=True starts everyone at the grid point closest
    to zero, mirroring the continuous lifecycle design.
    """
    rng = np.random.default_rng(seed)
    logy = np.empty((n_workers, n_years_keep))
    states, cums = {}, {}
    for name in ("1", "2"):
        grid, pi = disc[name]["grid"], disc[name]["pi"]
        if lifecycle:
            states[name] = np.full(n_workers, int(np.argmin(np.abs(grid))))
        else:
            states[name] = rng.choice(len(pi), size=n_workers, p=pi)
        cums[name] = np.cumsum(disc[name]["P"], axis=1)

    annual_flow = np.zeros(n_workers)
    year_idx = 0
    for q in range(4 * n_years_keep):
        for name in ("1", "2"):
            u = rng.random(n_workers)
            states[name] = (u[:, None] > cums[name][states[name]]).sum(axis=1)
        z = disc["1"]["grid"][states["1"]] + disc["2"]["grid"][states["2"]]
        annual_flow += np.exp(z)
        if q % 4 == 3:
            logy[:, year_idx] = np.log(annual_flow)
            annual_flow[:] = 0.0
            year_idx += 1
    return logy


def discretized_moments(params: dict, disc: dict | None = None, **sim_kwargs) -> dict:
    disc = disc or discretize_process(params)
    return moments_from_panel(simulate_discrete_panel(disc, **sim_kwargs))


def diagnose(params: dict, disc: dict) -> dict:
    """Per-component check that the chain's stationary variance matches theory."""
    rep = {}
    specs = {"1": (params["lambda1"], params["beta1"], params["sigma1"]),
             "2": (params["lambda2"], params["beta2"], params["sigma2"])}
    for name, (lam, beta, sig) in specs.items():
        g, pi = disc[name]["grid"], disc[name]["pi"]
        rep[name] = {"chain_var": float(pi @ (g - pi @ g) ** 2),
                     "theory_var": theoretical_stationary_var(lam, beta, sig),
                     "n_states": len(g)}
    c = disc["combined"]
    rep["combined"] = {"n_states": len(c["zgrid"]),
                       "rows_sum_to_1": bool(np.allclose(c["P"].sum(1), 1)),
                       "mean_earnings": float(c["pi"] @ c["egrid"])}
    return rep


def export_chain(disc: dict, outdir: str, prefix: str = "income_process") -> None:
    """
    Write the combined chain for step 2:
      *_grid.txt (earnings levels, mean 1), *_zgrid.txt (log grid),
      *_P.txt (quarterly transition matrix), *_pi.txt (stationary distribution).
    """
    import os
    c = disc["combined"]
    np.savetxt(os.path.join(outdir, f"{prefix}_grid.txt"), c["egrid"])
    np.savetxt(os.path.join(outdir, f"{prefix}_zgrid.txt"), c["zgrid"])
    np.savetxt(os.path.join(outdir, f"{prefix}_P.txt"), c["P"])
    np.savetxt(os.path.join(outdir, f"{prefix}_pi.txt"), c["pi"])
