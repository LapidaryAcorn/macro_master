"""
Sec 32 (v4): `chunk_workers` in simulate_log_annual_earnings.

1. chunk_workers=None (default) and chunk_workers >= n_workers are the
   unchunked panel, bit for bit.
2. Chunked panels are deterministic (common random numbers) and have the
   same moments as unchunked ones within Monte Carlo error (KMV design at the
   v4-type parameters, 8 seeds x 100k each side).
"""
from __future__ import annotations

import numpy as np

from kmv_earnings.simulate import MOMENT_ORDER, model_moments, simulate_log_annual_earnings
from kmv_earnings import kmv_method as km

P = {"lambda1": 0.113, "beta1": 0.99, "sigma1": 1.76,
     "lambda2": 0.0081, "beta2": 0.0092, "sigma2": 1.37}
SIM = dict(km.KMV_SIM, n_workers=20_000, exact_reset_burn=True)


def test_default_bit_identical():
    s = {k: v for k, v in SIM.items() if k not in ("cross_section", "d5_gap")}
    a = simulate_log_annual_earnings(P, **s)
    b = simulate_log_annual_earnings(P, chunk_workers=None, **s)
    c = simulate_log_annual_earnings(P, chunk_workers=20_000, **s)
    assert np.array_equal(a, b) and np.array_equal(a, c)


def test_chunked_deterministic_and_shape():
    s = {k: v for k, v in SIM.items() if k not in ("cross_section", "d5_gap")}
    a = simulate_log_annual_earnings(P, chunk_workers=6_000, **s)
    b = simulate_log_annual_earnings(P, chunk_workers=6_000, **s)
    assert a.shape == (20_000, 5) and np.array_equal(a, b)
    # blocks are distinct streams
    assert not np.array_equal(a[:6000], a[6000:12000])


def test_chunked_same_distribution():
    sim = dict(SIM, n_workers=100_000)
    un = np.array([[model_moments(P, seed=s, **sim)[k] for k in MOMENT_ORDER] for s in range(8)])
    ch = np.array([[model_moments(P, seed=100 + s, chunk_workers=25_000, **sim)[k]
                    for k in MOMENT_ORDER] for s in range(8)])
    se = np.sqrt(un.var(0, ddof=1) / 8 + ch.var(0, ddof=1) / 8)
    z = (ch.mean(0) - un.mean(0)) / se
    assert np.all(np.abs(z) < 3.5), dict(zip(MOMENT_ORDER, z.round(2)))
