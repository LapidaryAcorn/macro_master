"""
Sec 31: "the KMV method" options.

1. Defaults stay bit-for-bit identical to the code before Sec 31
   (simulate_log_annual_earnings, moments_from_panel, model_moments).
2. The new options do what KMV's code does: cross-section moments with a
   4-year "d5" reproduce KMV's own moments.txt from KMV's own simulated panel;
   kmv_generator reproduces KMV's published chain from their grid parameters.
   Both skip if the (local-only) replication package is absent.
3. dfls solves a least-squares problem.
"""
from __future__ import annotations

import os

import numpy as np
import pytest

from kmv_earnings.simulate import (KMV_TABLE3_PARAMS, MOMENT_ORDER, model_moments,
                                   moments_from_panel, simulate_log_annual_earnings)

R = os.path.join("data", "kmv_replication", "231376-V1 (2)", "20160042_data",
                 "Replication", "EarningsProcess")
HAVE_PKG = os.path.isdir(R)


def _simulate_pre_sec31(params, n_workers=50_000, n_years_keep=36, steps_per_quarter=6,
                        lifecycle=True, burn_in_years=0, seed=1234,
                        deterministic_arrival=False, reset_jumps=False):
    # verbatim copy of simulate_log_annual_earnings at commit 713f29f
    lam = np.array([params["lambda1"], params["lambda2"]]) / steps_per_quarter
    beta = np.array([params["beta1"], params["beta2"]]) / steps_per_quarter
    sig = np.array([params["sigma1"], params["sigma2"]])

    rng = np.random.default_rng(seed)
    decay = np.exp(-beta)
    p_jump = np.ones(2) if deterministic_arrival else 1.0 - np.exp(-lam)

    z = np.zeros((n_workers, 2))
    if not lifecycle:
        for j in range(2):
            if reset_jumps:
                sd = sig[j] * np.sqrt(lam[j] / (lam[j] + 2 * beta[j]))
            else:
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
            z[:, j] = np.where(u < p_jump[j], eps if reset_jumps else zd + eps, zd)
        if s >= burn_s:
            annual_flow += np.exp(z[:, 0] + z[:, 1])
            if (s - burn_s) % spy == spy - 1:
                logy[:, year_idx] = np.log(annual_flow)
                annual_flow[:] = 0.0
                year_idx += 1
    return logy


def _moments_pre_sec31(logy, n_age_cols=None):
    # verbatim copy of moments_from_panel at commit 713f29f
    if n_age_cols is not None:
        logy = logy[:, :n_age_cols]
    var_log = float(np.mean(np.var(logy, axis=0)))
    d1 = (logy[:, 1:] - logy[:, :-1]).ravel()
    d5 = (logy[:, 5:] - logy[:, :-5]).ravel()

    def kurt(x):
        x = x - x.mean()
        return float(np.mean(x**4) / np.mean(x**2) ** 2)

    return {
        "var_log_earns": var_log,
        "var_d1": float(np.var(d1)),
        "var_d5": float(np.var(d5)),
        "kurt_d1": kurt(d1),
        "kurt_d5": kurt(d5),
        "frac_d1_lt_10": float(np.mean(np.abs(d1) < 0.10)),
        "frac_d1_lt_20": float(np.mean(np.abs(d1) < 0.20)),
        "frac_d1_lt_50": float(np.mean(np.abs(d1) < 0.50)),
    }


P = KMV_TABLE3_PARAMS


@pytest.mark.parametrize("kw", [
    dict(n_workers=3_000, n_years_keep=12, seed=5),
    dict(n_workers=3_000, n_years_keep=12, seed=5, reset_jumps=True),
    dict(n_workers=3_000, n_years_keep=6, lifecycle=False, burn_in_years=2, seed=8),
    dict(n_workers=2_000, n_years_keep=6, seed=3, deterministic_arrival=True),
    dict(n_workers=2_000, n_years_keep=5, steps_per_quarter=4, lifecycle=False,
         burn_in_years=3, seed=4, reset_jumps=True),
])
def test_simulator_default_bit_identical(kw):
    a = simulate_log_annual_earnings(P, **kw)
    assert np.array_equal(a, _simulate_pre_sec31(P, **kw))
    assert np.array_equal(simulate_log_annual_earnings(P, kmv_euler=False, **kw), a)


def test_moments_default_bit_identical():
    y = simulate_log_annual_earnings(P, n_workers=4_000, n_years_keep=36, seed=2)
    for nac in (None, 31):
        old = _moments_pre_sec31(y, n_age_cols=nac)
        assert moments_from_panel(y, n_age_cols=nac) == old
        assert moments_from_panel(y, n_age_cols=nac, d5_gap=5, cross_section=False) == old
    kw = dict(n_workers=3_000, n_years_keep=36, seed=9)
    assert model_moments(P, n_age_cols=31, **kw) == _moments_pre_sec31(
        _simulate_pre_sec31(P, n_years_keep=31, n_workers=3_000, seed=9), n_age_cols=31)


def test_kmv_euler_differs():
    kw = dict(n_workers=3_000, n_years_keep=5, steps_per_quarter=4, lifecycle=False,
              burn_in_years=2, seed=4, reset_jumps=True)
    assert not np.array_equal(simulate_log_annual_earnings(P, kmv_euler=True, **kw),
                              simulate_log_annual_earnings(P, **kw))
    # no jumps (lambda = 0), start at 0 is not informative -> lifecycle with a
    # deterministic jump every step and no decay: Euler and exact coincide at beta=0
    q = dict(P, beta1=0.0, beta2=0.0)
    kw = dict(n_workers=1_000, n_years_keep=3, seed=1, deterministic_arrival=True)
    assert np.array_equal(simulate_log_annual_earnings(q, kmv_euler=True, **kw),
                          simulate_log_annual_earnings(q, **kw))


def test_cross_section_and_gap():
    rng = np.random.default_rng(0)
    y = rng.standard_normal((500, 6)).cumsum(1)
    m = moments_from_panel(y, cross_section=True, d5_gap=4)
    assert m["var_log_earns"] == pytest.approx(np.var(y[:, 0]))
    assert m["var_d5"] == pytest.approx(np.var(y[:, 4] - y[:, 0]))
    assert m["var_d1"] == pytest.approx(np.var(y[:, 1] - y[:, 0]))
    m4 = moments_from_panel(y, d5_gap=4)
    assert m4["var_d5"] == pytest.approx(np.var((y[:, 4:] - y[:, :-4]).ravel()))


@pytest.mark.skipif(not HAVE_PKG, reason="KMV replication package not present")
def test_moments_reproduce_kmv_output():
    E = os.path.join(R, "EarningsEstimation", "earnings_estimation_paper")
    y = np.loadtxt(os.path.join(E, "yannsim.txt"), delimiter=",")
    m = moments_from_panel(y, cross_section=True, d5_gap=4)
    kmv = {"var_log_earns": 0.695824263491323, "var_d1": 0.229240990711677,
           "var_d5": 0.464190169114165, "kurt_d1": 16.5302449008095,
           "kurt_d5": 12.0630368718262, "frac_d1_lt_10": 0.556999981403351,
           "frac_d1_lt_20": 0.672599971294403, "frac_d1_lt_50": 0.850600004196167}
    for k in MOMENT_ORDER:     # moments.txt prints ~15 digits; fractions are real(4)
        assert m[k] == pytest.approx(kmv[k], rel=1e-6), k


@pytest.mark.skipif(not HAVE_PKG, reason="KMV replication package not present")
def test_generator_reproduces_published_chain():
    from kmv_earnings.kmv_method import kmv_chain_from_grid_params
    from kmv_earnings.discretize import _combine
    D = os.path.join(R, "EarningsDiscretization", "earnings_discretization_paper")
    p = {"lambda1": 0.07987, "beta1": 0.7616, "sigma1": 1.7352,      # D/SetParameters l. 23-31
         "lambda2": 0.00656, "beta2": 0.00939, "sigma2": 1.5259}
    c = kmv_chain_from_grid_params(p, [3.78503497352302, 5.69800327154487],
                                   [0.792311684654520, 0.641191531992662])
    for k, n in (("1", 3), ("2", 11)):
        assert np.allclose(np.loadtxt(os.path.join(D, f"y{k}grid.txt")), c[k]["grid"], atol=1e-6)
        assert np.allclose(np.loadtxt(os.path.join(D, f"y{k}markov.txt")).reshape(n, n),
                           c[k]["Q"], atol=1e-6)
        assert np.allclose(np.loadtxt(os.path.join(D, f"y{k}dist.txt")), c[k]["pi"], atol=1e-6)
    disc = _combine(c)
    o = np.argsort(disc["combined"]["zgrid"], kind="stable")
    ym = np.loadtxt(os.path.join(D, "ymarkov_combined.txt")).reshape(33, 33)
    assert np.allclose(ym, disc["combined"]["Q"][np.ix_(o, o)], atol=1e-6)


def test_dfls_rosenbrock():
    from kmv_earnings.kmv_method import dfls
    res = lambda x: np.array([10 * (x[1] - x[0] ** 2), 1 - x[0]])
    x, f, nf = dfls(res, np.array([-1.2, 1.0]), rhobeg=0.5, rhoend=1e-8, maxfun=2000)
    assert f < 1e-10 and np.allclose(x, [1, 1], atol=1e-4)
