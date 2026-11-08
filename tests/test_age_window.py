"""
Sec 26 (A.2): the `n_age_cols` option of moments_from_panel / model_moments.
Default behaviour must stay bit-for-bit identical to the pre-Sec-26 code, and
the shortened simulation used by model_moments(n_age_cols=k) must give exactly
the moments of the full 36-year panel restricted to its first k columns.
"""
from __future__ import annotations

import numpy as np

from kmv_earnings.simulate import (KMV_TABLE3_PARAMS, moments_from_panel,
                                   model_moments, simulate_log_annual_earnings)


def _moments_pre_sec26(logy):
    # verbatim copy of moments_from_panel before the n_age_cols option
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


SIM = dict(n_workers=4_000, n_years_keep=36, seed=99)


def test_default_is_bit_identical():
    logy = simulate_log_annual_earnings(KMV_TABLE3_PARAMS, **SIM)
    assert moments_from_panel(logy) == _moments_pre_sec26(logy)
    assert model_moments(KMV_TABLE3_PARAMS, **SIM) == _moments_pre_sec26(logy)


def test_age_window_equals_column_slice():
    logy = simulate_log_annual_earnings(KMV_TABLE3_PARAMS, **SIM)
    ref = _moments_pre_sec26(logy[:, :31])
    assert moments_from_panel(logy, n_age_cols=31) == ref
    # model_moments simulates only 31 years; must equal the 36-year slice exactly
    assert model_moments(KMV_TABLE3_PARAMS, n_age_cols=31, **SIM) == ref


def test_age_window_changes_moments():
    logy = simulate_log_annual_earnings(KMV_TABLE3_PARAMS, **SIM)
    assert (moments_from_panel(logy, n_age_cols=31)["var_log_earns"]
            < moments_from_panel(logy)["var_log_earns"])
