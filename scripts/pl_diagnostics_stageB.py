"""
Stage B of the Poland diagnostics (see scripts/pl_diagnostics.py for Stage A
and output/pl_diagnostics/README.md). Needs a finished estimate in
output/<cc>/ for every country plotted -- run after Task 1's Poland
estimation (and, for the parameter-comparison panel, after Task 3's bootstrap
if CIs are wanted).

Produces, in output/pl_diagnostics/:
  05_model_fit_<cc>.png       -- 8 moments: data vs continuous vs discretized
  06_simulated_vs_grid_d1_<cc>.png -- simulated 1y-change ECDF vs GRID percentiles
  07_stationary_dist_<cc>.png -- chain stationary distribution, log and level
  08_sample_paths_<cc>.png    -- a few simulated 36-year worker paths
  09_var_by_age.png           -- var(log y) by age, all countries on one plot
  10_param_comparison.png     -- 6 params x 4 countries, KMV marked (+ bootstrap CIs if present)
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from kmv_earnings.grid_loader import _SCHEMAS, _detect_schema, _pct_prob_from_name
from kmv_earnings.simulate import (PARAM_ORDER, MOMENT_ORDER, MOMENT_LABELS,
                                    KMV_TABLE3_PARAMS, simulate_log_annual_earnings)
from kmv_earnings.discretize import discretize_process
from kmv_earnings.estimate import load_params, DEFAULT_BOUNDS

OUTDIR = os.path.join("output", "pl_diagnostics")
STYLE = {
    "USA": dict(color="#1b1b1b", marker="o", ls="-"),
    "FRA": dict(color="#4472c4", marker="s", ls="--"),
    "GER": dict(color="#c0392b", marker="^", ls=":"),
    "POL": dict(color="#1a9850", marker="D", ls="-."),
}
TARGETS_PATH = {cc: f"targets/{cc.lower()}_grid.json" for cc in STYLE}
DATA_FILES = {
    "USA": "data/Stats_20260720145608.csv", "FRA": "data/Stats_20260720145608.csv",
    "GER": "data/Stats_20260720145608.csv", "POL": "data/dane grid PL.csv",
}
plt.rcParams.update({"figure.dpi": 150, "savefig.dpi": 150, "font.size": 9,
                     "axes.grid": True, "grid.alpha": 0.3})


def available_countries():
    return [cc for cc in STYLE if os.path.exists(f"output/{cc.lower()}/params.json")]


def plot_model_fit(cc: str):
    tbl = pd.read_csv(f"output/{cc.lower()}/table2.csv")
    x = np.arange(len(tbl))
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    w = 0.25
    axes[0].bar(x - w, tbl["Data"], w, label="Data", color="#777")
    axes[0].bar(x, tbl["Model Estimated"], w, label="Model (continuous)", color="#4472c4")
    axes[0].bar(x + w, tbl["Model Discretized"], w, label="Model (discretized)", color="#c0392b")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(tbl["Moment"], rotation=45, ha="right", fontsize=7)
    axes[0].legend(fontsize=8)
    axes[0].set_title(f"{cc}: 8-moment fit")

    rel_c = (tbl["Model Estimated"] - tbl["Data"]) / tbl["Data"]
    rel_d = (tbl["Model Discretized"] - tbl["Data"]) / tbl["Data"]
    axes[1].axhline(0, color="k", lw=0.8)
    axes[1].bar(x - w / 2, rel_c * 100, w, label="continuous", color="#4472c4")
    axes[1].bar(x + w / 2, rel_d * 100, w, label="discretized", color="#c0392b")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(tbl["Moment"], rotation=45, ha="right", fontsize=7)
    axes[1].set_ylabel("% deviation from data")
    axes[1].legend(fontsize=8)
    axes[1].set_title(f"{cc}: relative deviation")
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, f"05_model_fit_{cc.lower()}.png"))
    plt.close(fig)


def plot_simulated_vs_grid_d1(cc: str):
    params = load_params(f"output/{cc.lower()}/params.json")
    logy = simulate_log_annual_earnings(params, n_workers=50_000, n_years_keep=36)
    d1 = (logy[:, 1:] - logy[:, :-1]).ravel()
    d1_sorted = np.sort(d1)
    ecdf = np.arange(1, len(d1_sorted) + 1) / len(d1_sorted)

    df = pd.read_csv(DATA_FILES[cc])
    df = df[(df["country"] == cc) & (df["gender"].astype(str).str.lower() == "male")]
    schema = _detect_schema(df)
    suffix = schema["d1_suffix"]
    pct_cols = [(c, p) for c in df.columns
                if (p := _pct_prob_from_name(c, suffix)) is not None
                and c not in (f"p9010{suffix}", f"p9050{suffix}", f"p5010{suffix}")]
    probs = np.array([p for _, p in pct_cols])
    order = np.argsort(probs)
    probs = probs[order]
    q = df[[pct_cols[i][0] for i in order]].mean(axis=0).to_numpy()

    fig, ax = plt.subplots(figsize=(7, 5))
    sub = (d1_sorted > np.percentile(d1_sorted, 0.1)) & (d1_sorted < np.percentile(d1_sorted, 99.9))
    ax.plot(d1_sorted[sub], ecdf[sub], color=STYLE[cc]["color"], lw=1.3,
            label="simulated ECDF (estimated process)")
    ax.scatter(q, probs, color="black", marker="x", s=35, zorder=5,
              label="GRID empirical percentiles (data)")
    ax.set_xlabel("1-year log change"); ax.set_ylabel("cumulative probability")
    ax.set_title(f"{cc}: simulated 1y-change distribution vs GRID percentiles")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, f"06_simulated_vs_grid_d1_{cc.lower()}.png"))
    plt.close(fig)


def plot_stationary_dist(cc: str):
    params = load_params(f"output/{cc.lower()}/params.json")
    disc = discretize_process(params)
    c = disc["combined"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    axes[0].stem(c["zgrid"], c["pi"], basefmt=" ", linefmt=STYLE[cc]["color"], markerfmt=".")
    axes[0].set_xlabel("log earnings grid point"); axes[0].set_ylabel("stationary probability")
    axes[0].set_title(f"{cc}: chain stationary distribution (log)")
    order = np.argsort(c["egrid"])
    axes[1].stem(c["egrid"][order], c["pi"][order], basefmt=" ", linefmt=STYLE[cc]["color"], markerfmt=".")
    axes[1].set_xlabel("earnings level (mean = 1)"); axes[1].set_ylabel("stationary probability")
    axes[1].set_title(f"{cc}: chain stationary distribution (level)")
    axes[1].set_xlim(0, np.percentile(c["egrid"], 99) * 1.2)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, f"07_stationary_dist_{cc.lower()}.png"))
    plt.close(fig)


def plot_sample_paths(cc: str, n_show: int = 6, seed: int = 55):
    params = load_params(f"output/{cc.lower()}/params.json")
    logy = simulate_log_annual_earnings(params, n_workers=2000, n_years_keep=36, seed=seed)
    rng = np.random.default_rng(seed)
    idx = rng.choice(logy.shape[0], size=n_show, replace=False)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    years = np.arange(1, 37)
    for i in idx:
        ax.plot(years, logy[i] - logy[i].mean(), lw=1.1, alpha=0.85)
    ax.axhline(0, color="grey", lw=0.6)
    ax.set_xlabel("year in panel (age proxy, years 1-36)")
    ax.set_ylabel("log annual earnings, demeaned per worker")
    ax.set_title(f"{cc}: {n_show} simulated 36-year worker paths\n"
                "(transitory wiggles + occasional persistent jumps)")
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, f"08_sample_paths_{cc.lower()}.png"))
    plt.close(fig)


def plot_var_by_age(countries: list[str]):
    fig, ax = plt.subplots(figsize=(7.5, 5))
    for cc in countries:
        params = load_params(f"output/{cc.lower()}/params.json")
        logy = simulate_log_annual_earnings(params, n_workers=50_000, n_years_keep=36)
        v = np.var(logy, axis=0)
        ax.plot(np.arange(1, 37), v, label=cc, lw=1.5, **{k: v2 for k, v2 in STYLE[cc].items() if k != "color"},
                color=STYLE[cc]["color"], ms=4)
        ax.axhline(v[-1], color=STYLE[cc]["color"], lw=0.6, ls=":", alpha=0.5)
    ax.set_xlabel("year in panel (proxy for labour-market age, 1-36)")
    ax.set_ylabel("var(log annual earnings)")
    ax.set_title("Dispersion build-up over the working life\n"
                "(dotted lines = year-36 level; distance from ergodic value not shown)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "09_var_by_age.png"))
    plt.close(fig)


def plot_param_comparison(countries: list[str]):
    rows = {cc: load_params(f"output/{cc.lower()}/params.json") for cc in countries}
    boot_path = "output/bootstrap/bootstrap_summary.csv"
    boot = pd.read_csv(boot_path) if os.path.exists(boot_path) else None

    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5))
    for ax, p in zip(axes.ravel(), PARAM_ORDER):
        xs = np.arange(len(countries))
        vals = [rows[cc][p] for cc in countries]
        colors = [STYLE[cc]["color"] for cc in countries]
        ax.bar(xs, vals, color=colors, width=0.6)
        if boot is not None:
            for i, cc in enumerate(countries):
                r = boot[boot.country == cc]
                if len(r) and pd.notna(r.iloc[0].get(f"{p}_se")):
                    se = r.iloc[0][f"{p}_se"]
                    ax.errorbar(i, vals[i], yerr=1.96 * se, color="black", capsize=4, lw=1.3)
        ax.axhline(KMV_TABLE3_PARAMS[p], color="grey", ls="--", lw=1.2, label="KMV (2018)")
        ax.set_xticks(xs); ax.set_xticklabels(countries)
        ax.set_title(p)
        if p == "lambda1":
            ax.legend(fontsize=7)
    fig.suptitle("Estimated parameters by country (KMV Table 3 marked)"
                + (" -- error bars = bootstrap 95% CI" if boot is not None else
                   " (no bootstrap CIs yet)"))
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(os.path.join(OUTDIR, "10_param_comparison.png"))
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--countries", nargs="+", default=None)
    args = ap.parse_args()
    countries = args.countries or available_countries()
    print(f"Stage B running for: {countries}")
    for cc in countries:
        plot_model_fit(cc)
        plot_simulated_vs_grid_d1(cc)
        plot_stationary_dist(cc)
        plot_sample_paths(cc)
    plot_var_by_age(countries)
    plot_param_comparison(countries)
    print(f"Stage B done -> {OUTDIR}")


if __name__ == "__main__":
    main()
