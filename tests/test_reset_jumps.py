"""
Sec 29.7: the `reset_jumps` option of simulate_log_annual_earnings.
Default (additive) must stay bit-for-bit identical to the pre-Sec-29.7 code,
for the lifecycle and the stationary-start design; reset must differ, start
from the reset stationary variance, and replace (not add to) the component.
"""
from __future__ import annotations

import numpy as np

from kmv_earnings.simulate import KMV_TABLE3_PARAMS, simulate_log_annual_earnings


def _simulate_pre_reset(params, n_workers=50_000, n_years_keep=36, steps_per_quarter=6,
                        lifecycle=True, burn_in_years=0, seed=1234,
                        deterministic_arrival=False):
    # verbatim copy of simulate_log_annual_earnings before the reset_jumps option
    lam = np.array([params["lambda1"], params["lambda2"]]) / steps_per_quarter
    beta = np.array([params["beta1"], params["beta2"]]) / steps_per_quarter
    sig = np.array([params["sigma1"], params["sigma2"]])

    rng = np.random.default_rng(seed)
    decay = np.exp(-beta)
    p_jump = np.ones(2) if deterministic_arrival else 1.0 - np.exp(-lam)

    z = np.zeros((n_workers, 2))
    if not lifecycle:
        for j in range(2):
            sd = sig[j] * np.sqrt((lam[j] / (2 * beta[j])) if beta[j] > 0 else 0.0)
            z[:, j] = rng.standard_normal(n_workers) * sd

    spy = 4 * steps_per_quarter
    total_s = spy * (burn_in_years + n_years_keep)
    burn_s = spy * burn_in_years

    logy = np.empty((n_workers, n_years_keep))
    annual_flow = np.zeros(n_workers)
    year_idx = 0

    for s in range(total_s):
        for j in range(2):
            u = rng.random(n_workers)
            eps = rng.standard_normal(n_workers) * sig[j]
            zd = z[:, j] * decay[j]
            z[:, j] = np.where(u < p_jump[j], zd + eps, zd)
        if s >= burn_s:
            annual_flow += np.exp(z[:, 0] + z[:, 1])
            if (s - burn_s) % spy == spy - 1:
                logy[:, year_idx] = np.log(annual_flow)
                annual_flow[:] = 0.0
                year_idx += 1
    return logy


P = KMV_TABLE3_PARAMS


def test_default_is_bit_identical_lifecycle():
    kw = dict(n_workers=3_000, n_years_keep=12, seed=5)
    assert np.array_equal(simulate_log_annual_earnings(P, **kw), _simulate_pre_reset(P, **kw))
    assert np.array_equal(simulate_log_annual_earnings(P, reset_jumps=False, **kw),
                          _simulate_pre_reset(P, **kw))


def test_default_is_bit_identical_stationary_and_ar1():
    kw = dict(n_workers=3_000, n_years_keep=6, lifecycle=False, burn_in_years=2, seed=8)
    assert np.array_equal(simulate_log_annual_earnings(P, **kw), _simulate_pre_reset(P, **kw))
    kw = dict(n_workers=2_000, n_years_keep=6, seed=3, deterministic_arrival=True)
    assert np.array_equal(simulate_log_annual_earnings(P, **kw), _simulate_pre_reset(P, **kw))


def test_reset_differs_and_replaces():
    kw = dict(n_workers=3_000, n_years_keep=8, seed=5)
    assert not np.array_equal(simulate_log_annual_earnings(P, reset_jumps=True, **kw),
                              simulate_log_annual_earnings(P, **kw))
    # with jumps every sub-step (p = 1) and no decay, reset leaves z = last eps:
    # annual log earnings then have the variance of a mean of iid lognormals,
    # whereas additive jumps accumulate into a random walk
    q = dict(P, beta1=1e-12, beta2=1e-12)
    yr = simulate_log_annual_earnings(q, n_workers=3_000, n_years_keep=6, seed=1,
                                      deterministic_arrival=True, reset_jumps=True)
    ya = simulate_log_annual_earnings(q, n_workers=3_000, n_years_keep=6, seed=1,
                                      deterministic_arrival=True)
    assert np.var(yr[:, -1]) < 0.5 * np.var(yr[:, 0]) + 1.0     # bounded
    assert np.var(ya[:, -1]) > 5 * np.var(ya[:, 0])              # growing


def test_reset_stationary_start_is_stationary():
    # lifecycle=False draws the start from N(0, lam sig^2/(lam + 2 beta)) under
    # reset; the cross-sectional variance must then stay flat over the years
    # (with the additive formula as the start it would fall toward it)
    y = simulate_log_annual_earnings(P, n_workers=100_000, n_years_keep=8, lifecycle=False,
                                     seed=2, reset_jumps=True)
    v = np.var(y, axis=0)
    assert abs(v[-1] / v[0] - 1) < 0.05


def test_reset_generator_lands_independently_of_origin():
    from kmv_earnings.discretize import component_grid_and_edges, component_generator
    g, e = component_grid_and_edges(11, 1.5)
    Qr = component_generator(g, e, 0.007, 0.0, 1.5, reset_jumps=True)   # beta=0: jumps only
    off = Qr - np.diag(np.diag(Qr))
    for j in range(11):            # rate into j is the same from every other origin
        col = np.delete(off[:, j], j)
        assert np.allclose(col, col[0])
    Qa = component_generator(g, e, 0.007, 0.0, 1.5)
    assert not np.allclose(Qa, Qr)
