"""
Sec 29.7: why did Sec 1 find that reset jumps miss KMV Table 2 (var d1 ~ 0.40)
while KMV's own code (reset) reproduces it?

Simulates KMV Table 3 parameters under several implementations and prints the
8 Table 2 moments for each:

  A  additive, lifecycle 36 yrs from z = 0         (our design, Sec 1-2)
  B  reset,    lifecycle 36 yrs from z = 0         (reset dropped into our design)
  C  reset,    stationary start + 25-yr burn, 5 yrs, true 5-yr change
  D  replica of KMV's EarningsEstimation code (description in DECISIONS 29.7):
     reset (rho = 0), Euler decay (1 - dt*beta), 16 steps/yr, start from
     N(0, lam sig^2/(2 beta + lam)), 25-yr burn, 5 annual obs, "5-yr" change =
     year 5 - year 1 (a 4-yr gap), one cross-section; KMV's unrounded estimates
  E  as D with Table 3 rounded parameters
  F  additive, stationary start + burn (ergodic), true 5-yr change
  H  as E but z frozen between jumps (no decay at all)

Usage: python scripts/reset_vs_sec1.py [--n 200000]
Output: output/pol_v3/kmv_method/reset_vs_sec1.json
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from kmv_earnings.simulate import (KMV_TABLE3_PARAMS, MOMENT_ORDER, moments_from_panel,
                                   simulate_log_annual_earnings)

# KMV's unrounded estimates (replication package,
# EarningsEstimation/earnings_estimation_paper/parameters.txt)
KMV_EXACT = {"lambda1": 0.07987252454850691, "beta1": 0.761606393760105,
             "sigma1": 1.73522438989472, "lambda2": 0.006560658609206301,
             "beta2": 0.009391785115626633, "sigma2": 1.52585125038250}
TABLE2_MODEL = {"var_log_earns": 0.70, "var_d1": 0.23, "var_d5": 0.46, "kurt_d1": 16.5,
                "kurt_d5": 12.1, "frac_d1_lt_10": 0.56, "frac_d1_lt_20": 0.67,
                "frac_d1_lt_50": 0.85}
# KMV's own printed model moments (earnings_estimation_paper/moments.txt)
KMV_CODE_OUTPUT = {"var_log_earns": 0.695824, "var_d1": 0.229241, "var_d5": 0.464190,
                   "kurt_d1": 16.5302, "kurt_d5": 12.0630, "frac_d1_lt_10": 0.5570,
                   "frac_d1_lt_20": 0.6726, "frac_d1_lt_50": 0.8506}


def kmv_replica(p, n, seed=7755, frozen=False):
    """Moments as KMV's Fortran computes them (one cross-section, 4-yr 'd5')."""
    rng = np.random.default_rng(seed)
    dt = 0.25                                   # quarters
    burn, keep = 100 * 4 + 1, 20 * 4 + 1        # steps (Tburn, Tsim - Tburn)
    lam = np.array([p["lambda1"], p["lambda2"]])
    beta = np.array([p["beta1"], p["beta2"]])
    sig = np.array([p["sigma1"], p["sigma2"]])
    z = rng.standard_normal((n, 2)) * np.sqrt(lam * sig**2 / (2 * beta + lam))
    levels = np.empty((n, keep))
    for t in range(1, burn + keep):
        jump = rng.random((n, 2)) < dt * lam
        eps = rng.standard_normal((n, 2)) * sig
        z = np.where(jump, eps, z if frozen else (1 - dt * beta) * z)
        if t >= burn:
            levels[:, t - burn] = np.exp(z.sum(1))
    ann = np.log(np.add.reduceat(levels[:, :80], np.arange(0, 80, 16), axis=1))  # 5 years
    d1, d5 = ann[:, 1] - ann[:, 0], ann[:, 4] - ann[:, 0]

    def k(x):
        x = x - x.mean(); return float(np.mean(x**4) / np.mean(x**2) ** 2)
    return {"var_log_earns": float(np.var(ann[:, 0])), "var_d1": float(np.var(d1)),
            "var_d5": float(np.var(d5)), "kurt_d1": k(d1), "kurt_d5": k(d5),
            "frac_d1_lt_10": float(np.mean(np.abs(d1) < .1)),
            "frac_d1_lt_20": float(np.mean(np.abs(d1) < .2)),
            "frac_d1_lt_50": float(np.mean(np.abs(d1) < .5))}


def kmv_objective(m, targets):
    w = {k: (0.5 if k.startswith("kurt") else 1.0) for k in MOMENT_ORDER}
    return sum(w[k] * (m[k] / targets[k] - 1) ** 2 for k in MOMENT_ORDER)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=200_000)
    a = ap.parse_args()
    T3, n = KMV_TABLE3_PARAMS, a.n
    life = dict(n_workers=n, n_years_keep=36, seed=11)
    erg = dict(n_workers=n, n_years_keep=5, lifecycle=False, burn_in_years=25, seed=11)
    res = {
        "A additive, lifecycle (ours)": moments_from_panel(simulate_log_annual_earnings(T3, **life)),
        "B reset, lifecycle from z=0": moments_from_panel(
            simulate_log_annual_earnings(T3, reset_jumps=True, **life)),
        "C reset, ergodic, 5-yr gap": None,
        "D KMV code replica, exact params": kmv_replica(KMV_EXACT, n),
        "E KMV code replica, Table 3": kmv_replica(T3, n),
        "F additive, ergodic, 5-yr gap": None,
        "H z frozen between jumps (no decay)": kmv_replica(T3, n, frozen=True),
    }
    # ergodic designs: one cross-section of 6 years so that a true 5-yr change exists
    for key, rj in (("C reset, ergodic, 5-yr gap", True), ("F additive, ergodic, 5-yr gap", False)):
        y = simulate_log_annual_earnings(T3, reset_jumps=rj, **dict(erg, n_years_keep=6))
        d1, d5 = y[:, 1] - y[:, 0], y[:, 5] - y[:, 0]

        def k(x):
            x = x - x.mean(); return float(np.mean(x**4) / np.mean(x**2) ** 2)
        res[key] = {"var_log_earns": float(np.var(y[:, 0])), "var_d1": float(np.var(d1)),
                    "var_d5": float(np.var(d5)), "kurt_d1": k(d1), "kurt_d5": k(d5),
                    "frac_d1_lt_10": float(np.mean(np.abs(d1) < .1)),
                    "frac_d1_lt_20": float(np.mean(np.abs(d1) < .2)),
                    "frac_d1_lt_50": float(np.mean(np.abs(d1) < .5))}
    res = {"KMV Table 2, model column": TABLE2_MODEL, "KMV code output (moments.txt)": KMV_CODE_OUTPUT,
           **res}
    data = {"var_log_earns": 0.70, "var_d1": 0.23, "var_d5": 0.46, "kurt_d1": 17.8,
            "kurt_d5": 11.55, "frac_d1_lt_10": 0.54, "frac_d1_lt_20": 0.71, "frac_d1_lt_50": 0.86}
    print(f"{'design':40s}" + "".join(f"{k[:10]:>11s}" for k in MOMENT_ORDER) + "   f_KMV")
    for name, m in res.items():
        print(f"{name:40s}" + "".join(f"{m[k]:11.3f}" for k in MOMENT_ORDER)
              + f"  {kmv_objective(m, data):7.4f}")
    os.makedirs("output/pol_v3/kmv_method", exist_ok=True)
    json.dump({"n": n, "kmv_exact_params": KMV_EXACT, "results": res},
              open("output/pol_v3/kmv_method/reset_vs_sec1.json", "w"), indent=2)


if __name__ == "__main__":
    main()
