"""
Sec 24, point 1: how much does measuring the level-dispersion target on a
3-year average (GRID's `std_perm`) instead of annual log earnings
(`std_loginc`, what moments_from_panel computes) move the MODEL moment?

Simulates the lifecycle panel at Poland's v1 canonical estimate
(output/pol/params.json), 50k workers, default estimation settings
(seed 1234, 36 years, 6 steps/quarter), and computes, per age (= panel
column), the cross-sectional variance of
  (a) annual log earnings                         -- moments_from_panel's var_log_earns
  (b) log of the 3-year average of earnings LEVELS over t-2..t
      -- GRID's permanent-earnings construction (codebook, std_perm);
         the model has no earnings threshold, so the ">= 2 qualifying
         years" rule never binds here.
averaged over ages. (b) can only start at the third panel year, so (a) is
also reported on that same column range for an apples-to-apples ratio.

Both are WITHIN-AGE variances (variance per column, then averaged), i.e.
already "age-residualized" in GRID's sense -- the model has no
deterministic age profile to remove.

Usage: python scripts/pol_perm3_check.py
Writes output/pol_v2/perm3_check.json.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from kmv_earnings.simulate import simulate_log_annual_earnings, moments_from_panel
from kmv_earnings.estimate import load_params


def var_log_perm3(logy: np.ndarray) -> float:
    y = np.exp(logy)
    perm = (y[:, 2:] + y[:, 1:-1] + y[:, :-2]) / 3.0      # years t-2..t, t = col 2..T-1
    return float(np.mean(np.var(np.log(perm), axis=0)))


def main():
    params = load_params("output/pol/params.json")
    logy = simulate_log_annual_earnings(params, n_workers=50_000, n_years_keep=36)
    out = {
        "params": params,
        "n_workers": 50_000,
        "var_log_annual_all_cols": moments_from_panel(logy)["var_log_earns"],
        "var_log_annual_cols_3plus": float(np.mean(np.var(logy[:, 2:], axis=0))),
        "var_log_perm3": var_log_perm3(logy),
    }
    out["ratio_perm3_to_annual_same_cols"] = (out["var_log_perm3"]
                                              / out["var_log_annual_cols_3plus"])
    for k, v in out.items():
        if k != "params":
            print(f"{k:34s} {v}")
    os.makedirs("output/pol_v2", exist_ok=True)
    with open("output/pol_v2/perm3_check.json", "w") as f:
        json.dump(out, f, indent=2)


if __name__ == "__main__":
    main()
