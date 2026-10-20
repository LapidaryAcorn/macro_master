"""
Sec 24: the two Polish diagnostics plots that depend on the level target or
the Polish estimate, regenerated for Poland v2. The v1 files are left
untouched; the new ones carry a _v2 suffix.

  01_moments_by_year_v2.png        -- as 01_moments_by_year.png, but POL's
                                      Var(log annual earnings) is GRID
                                      std_loginc^2 (year_panel_pol_v2.csv)
                                      instead of std_perm^2
  06_simulated_vs_grid_d1_pol_v2.png -- as 06_simulated_vs_grid_d1_pol.png,
                                      simulated at output/pol_v2/params.json

Plot code mirrors scripts/pl_diagnostics.py (stage A, figure 1) and
scripts/pl_diagnostics_stageB.py (plot_simulated_vs_grid_d1) exactly; only
the inputs and output names differ.

Usage: python scripts/pl_diagnostics_v2.py
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pl_diagnostics as A          # load_year_panel, STYLE, _shade_crisis, OUTDIR
import pl_diagnostics_stageB as B   # DATA_FILES, STYLE
from kmv_earnings.grid_loader import _detect_schema, _pct_prob_from_name
from kmv_earnings.simulate import simulate_log_annual_earnings
from kmv_earnings.estimate import load_params

OUTDIR = A.OUTDIR


def plot_moments_by_year_v2():
    panels = {cc: A.load_year_panel(cc) for cc in ("USA", "FRA", "GER")}
    panels["POL"] = pd.read_csv(os.path.join(OUTDIR, "year_panel_pol_v2.csv"))
    specs = [
        ("var_log", "Var(log annual earnings)"),
        ("var_d1", "Var(1y log change)"),
        ("var_d5", "Var(5y log change)"),
        ("kurt_d1", "Kurtosis(1y log change)"),
        ("kurt_d5", "Kurtosis(5y log change)"),
        ("nobs_d1", "N obs (1y-change cell)"),
    ]
    fig, axes = plt.subplots(3, 2, figsize=(10, 11), sharex=True)
    for ax, (key, title) in zip(axes.ravel(), specs):
        for cc, df in panels.items():
            ax.plot(df["year"], df[key], **A.STYLE[cc], ms=4, lw=1.3)
        A._shade_crisis(ax)
        ax.set_title(title)
        if key == "nobs_d1":
            ax.set_yscale("log")
            ax.yaxis.set_major_formatter(mticker.ScalarFormatter())
    axes[0, 0].legend(fontsize=7, loc="best")
    axes[-1, 0].set_xlabel("year")
    axes[-1, 1].set_xlabel("year")
    fig.suptitle("GRID moments by year, male prime-age (25-55): PL vs USA/FRA/GER\n"
                 "(grey band = 2008-10; PL var(log earnings) from std_loginc, Sec 24 v2)")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    path = os.path.join(OUTDIR, "01_moments_by_year_v2.png")
    fig.savefig(path)
    plt.close(fig)
    print("Saved:", path)


def plot_simulated_vs_grid_d1_pol_v2():
    cc = "POL"
    params = load_params("output/pol_v2/params.json")
    logy = simulate_log_annual_earnings(params, n_workers=50_000, n_years_keep=36)
    d1 = (logy[:, 1:] - logy[:, :-1]).ravel()
    d1_sorted = np.sort(d1)
    ecdf = np.arange(1, len(d1_sorted) + 1) / len(d1_sorted)

    df = pd.read_csv(B.DATA_FILES[cc])
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
    ax.plot(d1_sorted[sub], ecdf[sub], color=B.STYLE[cc]["color"], lw=1.3,
            label="simulated ECDF (estimated process, v2)")
    ax.scatter(q, probs, color="black", marker="x", s=35, zorder=5,
               label="GRID empirical percentiles (data)")
    ax.set_xlabel("1-year log change"); ax.set_ylabel("cumulative probability")
    ax.set_title(f"{cc} v2: simulated 1y-change distribution vs GRID percentiles")
    ax.legend(fontsize=8)
    fig.tight_layout()
    path = os.path.join(OUTDIR, "06_simulated_vs_grid_d1_pol_v2.png")
    fig.savefig(path)
    plt.close(fig)
    print("Saved:", path)


if __name__ == "__main__":
    plot_moments_by_year_v2()
    plot_simulated_vs_grid_d1_pol_v2()
