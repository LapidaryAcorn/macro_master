"""DECISIONS Sec 33: estimate_dfls(upper=...) -- the default box is unchanged,
and a wider beta box (upper 4, v4_nob) maps and polishes correctly."""
import numpy as np

from kmv_earnings import kmv_method as km

SIM = dict(km.KMV_SIM, n_workers=3000, seed=11)
GUESS = {"lambda1": 0.15, "beta1": 0.8, "sigma1": 1.7,
         "lambda2": 0.01, "beta2": 0.008, "sigma2": 1.2}
UP4 = dict(km.KMV_UPPER, beta1=4.0, beta2=4.0)


def _run(upper):
    return km.estimate_dfls(km.KMV_TARGETS, sim_kwargs=SIM, guess=GUESS, ndfls=1,
                            rhobeg=1.0, maxfun=25, upper=upper)


def test_default_box_bit_identical():
    p0, i0 = _run(None)
    p1, i1 = _run(dict(km.KMV_UPPER))
    assert p0 == p1
    assert i0["f"] == i1["f"] and i0["n_evals"] == i1["n_evals"]


def test_wide_box_roundtrip_and_beyond_one():
    p = dict(GUESS, beta1=2.5, beta2=0.3)
    x = km._to_x(p, upper=UP4)
    q = km._to_p(x, upper=UP4)
    assert all(np.isclose(p[k], q[k], rtol=1e-12) for k in p)
    # same parameter vector, different box -> different logistic coordinates
    assert not np.allclose(km._to_x(GUESS), km._to_x(GUESS, upper=UP4))


def test_wide_box_polish_runs_and_respects_box():
    p, info = _run(UP4)
    assert 0 < p["beta1"] < 4.0 and 0 < p["beta2"] < 4.0
    assert np.isfinite(info["f"]) and info["n_evals"] <= 25
