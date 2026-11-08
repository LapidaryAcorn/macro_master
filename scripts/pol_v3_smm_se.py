"""
Sec 27 (C): asymptotic SMM standard errors at the canonical Poland v3 estimate,
as a check on the bootstrap (Sec 26, B.4).

    V = (G'WG)^-1 G'W (Omega_data + Omega_sim) W G (G'WG)^-1

  G           d m / d theta (8 x 6), central differences in log theta on a
              large CRN panel (200k workers, seed 1234, ages 25-55), at two
              step sizes (h = 0.02 and 0.05 in logs) to check stability;
              converted to d m / d theta by dividing by theta.
  W           diag(1 / m_hat^2), m_hat = v3 targets -- the weights of our
              objective (sum of squared relative deviations).
  Omega_data  covariance of the annual moments around their linear trends
              (output/pol_v3/trends/trend_residuals.csv) divided by the number
              of years, because each target is a mean over years. Series have
              different lengths (d5: 16 yrs, d1: 20, level: 21), so element
              (j,k) = s_jk * n_jk / (T_j T_k), s_jk = sum over the n_jk common
              years of r_j r_k / (n_jk - 2) (OLS residuals, 2 parameters).
              Also the non-detrended version (deviations from the mean,
              divisor n_jk - 1). This is year-to-year variability of the GRID
              statistics, NOT sampling error from microdata.
  Omega_sim   covariance of the model moments across 60 simulation seeds on
              the estimation-size panel (50k workers, ages 25-55).

Derived quantities (var2, ergodic sd of z2, half-life of z2) by the delta
method with numerical gradients.

Usage: python scripts/pol_v3_smm_se.py [--params output/pol_v3/params.json] [--workers 10]
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from kmv_earnings.simulate import MOMENT_ORDER, PARAM_ORDER, model_moments

N_AGE = 31
TARGETS = "targets/pol_grid_v3.json"
RESID = "output/pol_v3/trends/trend_residuals.csv"
PANEL = "output/pol_v3/year_panel_pol_v3.csv"
OUT = "output/pol_v3/smm_se"
G_SIM = dict(n_workers=int(os.environ.get("SMM_SE_GN", 200_000)), n_years_keep=36, n_age_cols=N_AGE, seed=1234)
SIM_SEEDS = list(range(5000, 5000 + int(os.environ.get("SMM_SE_NSEEDS", 60))))
COLS = {"var_log_earns": "var_log", "var_d1": "var_d1", "var_d5": "var_d5",
        "kurt_d1": "kurt_d1", "kurt_d5": "kurt_d5", "frac_d1_lt_10": "frac_d1_lt_10",
        "frac_d1_lt_20": "frac_d1_lt_20", "frac_d1_lt_50": "frac_d1_lt_50"}


def _mvec(args):
    p, sim = args
    m = model_moments(p, **sim)
    return np.array([m[k] for k in MOMENT_ORDER])


def jacobian(p, h, pool):
    jobs = []
    for k in PARAM_ORDER:
        for s in (+1, -1):
            q = dict(p); q[k] = p[k] * np.exp(s * h)
            jobs.append((q, G_SIM))
    out = pool.map(_mvec, jobs)
    G = np.empty((8, 6))
    for i, k in enumerate(PARAM_ORDER):
        G[:, i] = (out[2 * i] - out[2 * i + 1]) / (2 * h) / p[k]     # d m / d theta
    return G


def omega_data(df: pd.DataFrame, ddof: int) -> np.ndarray:
    R = df[MOMENT_ORDER].to_numpy()
    T = np.sum(~np.isnan(R), axis=0)
    O = np.empty((8, 8))
    for j in range(8):
        for k in range(8):
            both = ~np.isnan(R[:, j]) & ~np.isnan(R[:, k])
            n = both.sum()
            s = np.sum(R[both, j] * R[both, k]) / (n - ddof)
            O[j, k] = s * n / (T[j] * T[k])
    return O


def sandwich(G, W, Om):
    A = np.linalg.inv(G.T @ W @ G)
    return A @ G.T @ W @ Om @ W @ G @ A


def derived(theta):
    p = dict(zip(PARAM_ORDER, theta))
    v2 = p["lambda2"] * p["sigma2"] ** 2 / (2 * p["beta2"])
    return np.array([v2, np.sqrt(v2), np.log(2) / p["beta2"] / 4])


DERIVED = ["var2", "erg_sd_z2", "half_life_2_yr"]


def delta(theta, V):
    g0 = derived(theta)
    J = np.empty((len(g0), 6))
    for i in range(6):
        e = np.zeros(6); e[i] = theta[i] * 1e-5
        J[:, i] = (derived(theta + e) - derived(theta - e)) / (2 * e[i])
    return g0, np.sqrt(np.diag(J @ V @ J.T))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", default="output/pol_v3/params.json")
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    p = json.load(open(a.params))
    theta = np.array([p[k] for k in PARAM_ORDER])
    t = json.load(open(TARGETS))
    mhat = np.array([t[k] for k in MOMENT_ORDER])
    W = np.diag(1 / mhat ** 2)

    with Pool(a.workers) as pool:
        Gs = {h: jacobian(p, h, pool) for h in (0.02, 0.05)}
        sims = np.array(pool.map(_mvec, [(p, dict(n_workers=50_000, n_years_keep=36,
                                                   n_age_cols=N_AGE, seed=s)) for s in SIM_SEEDS]))
    Om_sim = np.cov(sims.T)
    res = pd.read_csv(RESID)
    panel = pd.read_csv(PANEL)
    dev = panel[["year"]].copy()
    for k, c in COLS.items():
        dev[k] = panel[c] - panel[c].mean()
    Om_dt = omega_data(res, ddof=2)
    Om_raw = omega_data(dev, ddof=1)

    G = Gs[0.02]
    rel_G = np.abs(Gs[0.05] - Gs[0.02]) / np.maximum(np.abs(Gs[0.02]), 1e-12)
    cases = {"detrended+sim": Om_dt + Om_sim, "detrended": Om_dt,
             "raw+sim": Om_raw + Om_sim, "raw": Om_raw, "sim_only": Om_sim}
    rows = []
    for name, Om in cases.items():
        for h, GG in Gs.items():
            V = sandwich(GG, W, Om)
            se = np.sqrt(np.diag(V))
            g0, gse = delta(theta, V)
            row = {"omega": name, "h": h}
            row.update({f"se_{k}": se[i] for i, k in enumerate(PARAM_ORDER)})
            row.update({f"relse_{k}": se[i] / theta[i] for i, k in enumerate(PARAM_ORDER)})
            row.update({f"se_{k}": gse[i] for i, k in enumerate(DERIVED)})
            rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(a.out, "smm_se.csv"), index=False)
    g0, _ = delta(theta, np.zeros((6, 6)))
    out = {"params": p, "derived": dict(zip(DERIVED, g0)), "G_h0.02": Gs[0.02].tolist(),
           "G_h0.05": Gs[0.05].tolist(), "G_max_rel_diff_by_param":
           dict(zip(PARAM_ORDER, rel_G.max(0).tolist())),
           "G_median_rel_diff": float(np.median(rel_G)),
           "Omega_sim": Om_sim.tolist(), "Omega_data_detrended": Om_dt.tolist(),
           "Omega_data_raw": Om_raw.tolist(),
           "sd_moments": {"data_detrended": dict(zip(MOMENT_ORDER, np.sqrt(np.diag(Om_dt)))),
                          "data_raw": dict(zip(MOMENT_ORDER, np.sqrt(np.diag(Om_raw)))),
                          "sim": dict(zip(MOMENT_ORDER, np.sqrt(np.diag(Om_sim))))},
           "cond_GWG": float(np.linalg.cond(G.T @ W @ G)),
           "rows": rows}
    json.dump(out, open(os.path.join(a.out, "smm_se.json"), "w"), indent=2, default=float)

    pd.set_option("display.width", 250)
    print("sd of each moment (data detrended / raw / sim):")
    for k in MOMENT_ORDER:
        print(f"  {k:15s} {out['sd_moments']['data_detrended'][k]:.5g}"
              f"  {out['sd_moments']['data_raw'][k]:.5g}  {out['sd_moments']['sim'][k]:.5g}")
    print("G stability, max rel. diff h=0.05 vs 0.02 by parameter column:",
          {k: round(v, 3) for k, v in out["G_max_rel_diff_by_param"].items()},
          " median", round(out["G_median_rel_diff"], 4))
    print(df[[c for c in df.columns if c.startswith(("omega", "h", "relse_", "se_var2",
                                                     "se_erg", "se_half"))]]
          .to_string(index=False, float_format="%.4g"))


if __name__ == "__main__":
    main()
