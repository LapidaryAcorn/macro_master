"""Exact time aggregation of the reset jump-drift process (DECISIONS Sec 33).

Continuous time, two components j = 1, 2 (rates per quarter):
  Poisson arrivals at rate lambda_j; on arrival z_j is RESET to eps ~ N(0, sigma_j^2);
  between arrivals z_j(t) = z_j(tau) exp(-beta_j (t - tau)).
Annual income = integral over the year of exp(z1(t) + z2(t)) dt; the moments
use its log (the constant time unit drops out).

Integration over a step of length h quarters:
  * z1 exactly: on a segment of length T with start value a,
        int_0^T exp(a e^{-beta t}) dt = [Ei(a) - Ei(a e^{-beta T})] / beta
    (substitute u = a e^{-beta t}); Ei = scipy.special.expi.
  * z2 is held at its value within the step (its decay is applied at the end
    of the step), but its jumps -- and those of z1 -- are placed at their
    exact times inside the step, so the year is split exactly at every jump.
  * Up to `max_jumps` arrivals per component and step (Poisson count by
    inversion, capped); with lambda*h ~ 0.04 the cap is never binding in practice.
Common random numbers: every worker draws the same fixed block of uniforms
and normals at every step whatever the parameters (count by inversion, slot
times and slot shocks), so the panel is a deterministic, piecewise-smooth
function of the parameters, as in the Euler simulator.

Start: the stationary normal, then the exact law of the state after the
burn-in (as `exact_reset_burn`): z = eps e^{-beta K} if a jump occurred in the
burn-in, K = time since the last one (truncated exponential), else
z0 e^{-beta B}.
"""
from __future__ import annotations

import numpy as np
from scipy.special import expi

_GAMMA = 0.5772156649015329
_SMALL = 1e-8


def seg_integral(a, beta: float, T):
    """int_0^T exp(a exp(-beta t)) dt, vectorized in a and T (T >= 0)."""
    a = np.asarray(a, dtype=float)
    T = np.asarray(T, dtype=float)
    bt = beta * T
    ae = a * np.exp(-bt)
    out = np.empty(np.broadcast(a, T).shape)
    small_a = np.abs(a) < _SMALL
    small_e = (np.abs(ae) < _SMALL) & ~small_a
    big = ~(small_a | small_e)
    a_b, ae_b = np.broadcast_to(a, out.shape), np.broadcast_to(ae, out.shape)
    T_b, bt_b = np.broadcast_to(T, out.shape), np.broadcast_to(bt, out.shape)
    if big.any():
        out[big] = (expi(a_b[big]) - expi(ae_b[big])) / beta
    if small_a.any():                        # exp(.) ~ 1 + a e^{-beta t}
        out[small_a] = T_b[small_a] + (a_b[small_a] - ae_b[small_a]) / beta
    if small_e.any():                        # Ei(x) ~ gamma + ln|x| + x, ln|ae| = ln|a| - beta T
        x = a_b[small_e]
        ei_e = _GAMMA + np.log(np.abs(x)) - bt_b[small_e] + ae_b[small_e]
        out[small_e] = (expi(x) - ei_e) / beta
    return out


def _poisson_count(u, mu: float, kmax: int):
    """Poisson(mu) count by inversion of u, capped at kmax."""
    n = np.zeros(u.shape, dtype=np.int64)
    p = np.exp(-mu)
    c = p
    for k in range(1, kmax + 1):
        n += u > c
        p *= mu / k
        c += p
    return n


def simulate_exact(params: dict, n_workers: int = 50_000, n_years_keep: int = 5,
                   burn_in_years: int = 25, seed=1234, steps_per_quarter: int = 1,
                   max_jumps: int = 3, chunk_workers: int | None = None, **_ignored) -> np.ndarray:
    """Log annual earnings (n_workers, n_years_keep) of the continuous reset
    process with exact time aggregation; params = continuous-time rates per
    quarter (lambda_c, beta_c). steps_per_quarter = M sets the step h = 1/M
    quarter on which z2 is held fixed (z1 and all jump times are exact)."""
    if chunk_workers is not None and n_workers > chunk_workers:
        sizes = [chunk_workers] * (n_workers // chunk_workers)
        if n_workers % chunk_workers:
            sizes.append(n_workers % chunk_workers)
        kw = dict(n_years_keep=n_years_keep, burn_in_years=burn_in_years,
                  steps_per_quarter=steps_per_quarter, max_jumps=max_jumps)
        return np.vstack([simulate_exact(params, n_workers=n, seed=[seed, i], **kw)
                          for i, n in enumerate(sizes)])
    lam = np.array([params["lambda1"], params["lambda2"]], dtype=float)
    beta = np.array([params["beta1"], params["beta2"]], dtype=float)
    sig = np.array([params["sigma1"], params["sigma2"]], dtype=float)
    W, K, M = n_workers, max_jumps, steps_per_quarter
    h = 1.0 / M
    rng = np.random.default_rng(seed)

    # stationary normal start, then the exact law after the burn-in
    z = np.empty((W, 2))
    for j in range(2):
        z[:, j] = rng.standard_normal(W) * sig[j] * np.sqrt(lam[j] / (lam[j] + 2 * beta[j]))
    B = 4.0 * burn_in_years
    if B > 0:
        for j in range(2):
            u = rng.random(W)
            eps = rng.standard_normal(W) * sig[j]
            q = np.exp(-lam[j] * B)
            jumped = u >= q
            v = np.where(jumped, (u - q) / (1.0 - q), 0.0)
            age = -np.log1p(-v * (1.0 - q)) / lam[j]
            z[:, j] = np.where(jumped, eps * np.exp(-beta[j] * age), z[:, j] * np.exp(-beta[j] * B))
    a, b = z[:, 0].copy(), z[:, 1].copy()
    d1, d2 = np.exp(-beta[0] * h), np.exp(-beta[1] * h)

    spy = 4 * M
    logy = np.empty((W, n_years_keep))
    flow = np.zeros(W)
    for s in range(spy * n_years_keep):
        U = rng.random((W, 2, 1 + K))
        E = rng.standard_normal((W, 2, K))
        n1 = _poisson_count(U[:, 0, 0], lam[0] * h, K)
        n2 = _poisson_count(U[:, 1, 0], lam[1] * h, K)
        ev = (n1 + n2) > 0
        quiet = ~ev
        # no arrival in the step: one exact segment
        aq = a[quiet]
        flow[quiet] += np.exp(b[quiet]) * seg_integral(aq, beta[0], h)
        a[quiet] = aq * d1
        if ev.any():
            idx = np.flatnonzero(ev)
            slot = np.arange(K)
            t1 = np.where(slot < n1[idx, None], h * U[idx, 0, 1:], h)      # unused slots at the end
            t2 = np.where(slot < n2[idx, None], h * U[idx, 1, 1:], h)
            tt = np.concatenate([t1, t2], axis=1)
            typ = np.concatenate([np.where(slot < n1[idx, None], 1, 0),
                                  np.where(slot < n2[idx, None], 2, 0)], axis=1)
            ee = np.concatenate([E[idx, 0, :] * sig[0], E[idx, 1, :] * sig[1]], axis=1)
            o = np.argsort(tt, axis=1, kind="stable")
            tt = np.take_along_axis(tt, o, 1)
            typ = np.take_along_axis(typ, o, 1)
            ee = np.take_along_axis(ee, o, 1)
            aS, bS = a[idx], b[idx]
            fS = np.zeros(idx.size)
            s0 = np.zeros(idx.size)
            for k in range(2 * K):
                dt = tt[:, k] - s0
                fS += np.exp(bS) * seg_integral(aS, beta[0], dt)
                aS = aS * np.exp(-beta[0] * dt)
                s0 = tt[:, k]
                aS = np.where(typ[:, k] == 1, ee[:, k], aS)
                bS = np.where(typ[:, k] == 2, ee[:, k], bS)
            dt = h - s0
            fS += np.exp(bS) * seg_integral(aS, beta[0], dt)
            aS = aS * np.exp(-beta[0] * dt)
            flow[idx] += fS
            a[idx], b[idx] = aS, bS
        b *= d2
        if s % spy == spy - 1:
            logy[:, s // spy] = np.log(flow)
            flow[:] = 0.0
    return logy
