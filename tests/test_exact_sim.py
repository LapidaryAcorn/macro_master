"""DECISIONS Sec 33: exact time aggregation (kmv_earnings/exact_sim.py)."""
import numpy as np
from scipy.integrate import quad

from kmv_earnings.exact_sim import seg_integral, simulate_exact
from kmv_earnings.simulate import simulate_log_annual_earnings

P = {"lambda1": 0.164, "beta1": 3.07, "sigma1": 1.82,
     "lambda2": 0.0123, "beta2": 0.0065, "sigma2": 1.12}


def test_segment_integral_matches_quadrature():
    for a in (-6.0, -1.0, -1e-9, 0.0, 1e-9, 0.3, 4.0, 8.0):
        for beta in (0.006, 1.15, 3.07, 20.0):
            for T in (0.0, 0.01, 0.25, 1.0):
                ex = quad(lambda t: np.exp(a * np.exp(-beta * t)), 0, T,
                          epsabs=1e-14, epsrel=1e-12)[0] if T > 0 else 0.0
                got = seg_integral(np.array([a]), beta, np.array([T]))[0]
                assert abs(got - ex) <= 1e-9 * max(1.0, abs(ex)), (a, beta, T)


def test_mean_flow_matches_analytic_z1_only():
    q = dict(P, sigma2=1e-6)
    lam, b, s = q["lambda1"], q["beta1"], q["sigma1"]
    ez = quad(lambda A: lam * np.exp(-lam * A) * np.exp(s * s * np.exp(-2 * b * A) / 2),
              0, np.inf, limit=200)[0]
    y = simulate_exact(q, n_workers=100_000, n_years_keep=3, seed=5)
    assert abs(np.exp(y).mean() / (4 * ez) - 1) < 0.01


def test_deterministic_and_routed():
    kw = dict(n_workers=2000, n_years_keep=5, burn_in_years=25, seed=3)
    y1 = simulate_exact(P, **kw)
    y2 = simulate_log_annual_earnings(P, lifecycle=False, reset_jumps=True,
                                      exact_integral=True, steps_per_quarter=1, **kw)
    assert np.array_equal(y1, y2) and np.isfinite(y1).all()
