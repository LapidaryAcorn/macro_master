"""
Diagnostic plots for the Polish GRID data and the estimated earnings process,
to look at before trusting the estimation (requested alongside Task 1,
2026-10-11). Two stages:

  Stage A (raw data only, no estimation needed):
    - year-by-year panels of all 8 moments + nobs, PL vs USA/FRA/GER
    - quantile-vs-normal shape plots for selected years

  Stage B (needs a finished Poland estimate in output/pol/):
    - 8-moment fit: data vs continuous model vs discretized chain
    - simulated 1y-change distribution vs empirical GRID percentiles
    - chain stationary distribution, sample worker paths, var(log y) by age
    - 6-parameter cross-country comparison (+ bootstrap CIs if available)

Usage:
    python scripts/pl_diagnostics.py --stage A
    python scripts/pl_diagnostics.py --stage B
    python scripts/pl_diagnostics.py --stage all
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
import matplotlib.ticker as mticker

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from kmv_earnings.grid_loader import (_SCHEMAS, _detect_schema,
                                       _pct_prob_from_name,
                                       fractions_from_percentiles)
from kmv_earnings.simulate import (PARAM_ORDER, KMV_TABLE3_PARAMS,
                                    simulate_log_annual_earnings)
from kmv_earnings.discretize import (discretize_process, theoretical_stationary_var)
from kmv_earnings.estimate import load_params

OUTDIR = os.path.join("output", "pl_diagnostics")
os.makedirs(OUTDIR, exist_ok=True)

DATA_FILES = {
    "USA": "data/Stats_20260720145608.csv",
    "FRA": "data/Stats_20260720145608.csv",
    "GER": "data/Stats_20260720145608.csv",
    "POL": "data/dane grid PL.csv",
}
LEVEL_NOBS = {"res_1yr_log_chg": "nobs_log_inc", "resid_dNyr": "nobs_perm"}

# Consistent per-country style: distinct colour AND marker/linestyle so the
# figures remain legible printed in black and white.
STYLE = {
    "USA": dict(color="#1b1b1b", marker="o", ls="-",  label="USA"),
    "FRA": dict(color="#4472c4", marker="s", ls="--", label="FRA"),
    "GER": dict(color="#c0392b", marker="^", ls=":",  label="GER"),
    "POL": dict(color="#1a9850", marker="D", ls="-.", label="POL"),
}
CRISIS = (2008, 2010)

plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 150, "font.size": 9,
    "axes.grid": True, "grid.alpha": 0.3,
})


def _shade_crisis(ax):
    ax.axvspan(CRISIS[0], CRISIS[1], color="grey", alpha=0.15, zorder=0,
               label="2008-10 crisis" if ax.get_legend() is None else None)


def load_year_panel(cc: str) -> pd.DataFrame:
    """Per-year moments + nobs for one country, NOT averaged over years."""
    df = pd.read_csv(DATA_FILES[cc])
    df = df[df["country"] == cc]
    df = df[df["gender"].astype(str).str.lower() == "male"]
    assert df["gender"].nunique() == 1
    schema = _detect_schema(df)
    rows = []
    suffix = schema["d1_suffix"]
    pct_cols = [(c, p) for c in df.columns
                if (p := _pct_prob_from_name(c, suffix)) is not None
                and c not in (f"p9010{suffix}", f"p9050{suffix}", f"p5010{suffix}")]
    probs = np.array([p for _, p in pct_cols])
    order = np.argsort(probs)
    probs = probs[order]
    pct_col_names = [pct_cols[i][0] for i in order]
    for _, r in df.sort_values("year").iterrows():
        row = {
            "year": int(r["year"]),
            "var_log": r[schema["level_std"]] ** 2 if pd.notna(r[schema["level_std"]]) else np.nan,
            "var_d1": r[schema["d1_std"]] ** 2 if pd.notna(r[schema["d1_std"]]) else np.nan,
            "var_d5": r[schema["d5_std"]] ** 2 if pd.notna(r[schema["d5_std"]]) else np.nan,
            "kurt_d1": r[schema["d1_kurt"]],
            "kurt_d5": r[schema["d5_kurt"]],
            "nobs_level": r.get(LEVEL_NOBS[schema["name"]], np.nan),
            "nobs_d1": r[schema["d1_nobs"]],
            "nobs_d5": r[schema["d5_nobs"]],
        }
        q = r[pct_col_names].to_numpy(dtype=float)
        if np.isfinite(q).sum() >= 5:
            fracs = fractions_from_percentiles(probs[np.isfinite(q)], q[np.isfinite(q)])
        else:
            fracs = {k: np.nan for k in ("frac_d1_lt_10", "frac_d1_lt_20", "frac_d1_lt_50")}
        row.update(fracs)
        rows.append(row)
    return pd.DataFrame(rows)


def stage_a():
    panels = {cc: load_year_panel(cc) for cc in ("USA", "FRA", "GER", "POL")}
    for cc, df in panels.items():
        df.to_csv(os.path.join(OUTDIR, f"year_panel_{cc.lower()}.csv"), index=False)

    # ---- 1. main time-series grid ----
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
            ax.plot(df["year"], df[key], **STYLE[cc], ms=4, lw=1.3)
        _shade_crisis(ax)
        ax.set_title(title)
        if key == "nobs_d1":
            ax.set_yscale("log")
            ax.yaxis.set_major_formatter(mticker.ScalarFormatter())
    axes[0, 0].legend(fontsize=7, loc="best")
    axes[-1, 0].set_xlabel("year")
    axes[-1, 1].set_xlabel("year")
    fig.suptitle("GRID moments by year, male prime-age (25-55): PL vs USA/FRA/GER\n"
                 "(grey band = 2008-10)")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(os.path.join(OUTDIR, "01_moments_by_year.png"))
    plt.close(fig)

    # ---- 2. fractions of small changes ----
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharex=True, sharey=True)
    for ax, key, title in zip(
            axes, ("frac_d1_lt_10", "frac_d1_lt_20", "frac_d1_lt_50"),
            ("P(|d1y log chg| < 10%)", "< 20%", "< 50%")):
        for cc, df in panels.items():
            ax.plot(df["year"], df[key], **STYLE[cc], ms=4, lw=1.3)
        _shade_crisis(ax)
        ax.set_title(title)
        ax.set_xlabel("year")
    axes[0].set_ylabel("fraction of workers")
    axes[0].legend(fontsize=8)
    fig.suptitle("Fractions of small 1-year changes, by year (grey band = 2008-10)")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(os.path.join(OUTDIR, "02_fractions_by_year.png"))
    plt.close(fig)

    # ---- 3. nobs / level sample size ----
    fig, ax = plt.subplots(figsize=(7, 4.2))
    for cc, df in panels.items():
        ax.plot(df["year"], df["nobs_level"], **STYLE[cc], ms=4, lw=1.3)
    _shade_crisis(ax)
    ax.set_yscale("log")
    ax.yaxis.set_major_formatter(mticker.ScalarFormatter())
    ax.set_title("Sample size (level-earnings cell), by year")
    ax.set_xlabel("year"); ax.set_ylabel("N obs (log scale)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "03_nobs_level_by_year.png"))
    plt.close(fig)

    # ---- 4. quantile-vs-normal shape plots ----
    def qq_panel(ax, cc, year, color):
        df = pd.read_csv(DATA_FILES[cc])
        df = df[(df["country"] == cc) & (df["year"] == year)]
        df = df[df["gender"].astype(str).str.lower() == "male"]
        schema = _detect_schema(df)
        suffix = schema["d1_suffix"]
        pct_cols = [(c, p) for c in df.columns
                    if (p := _pct_prob_from_name(c, suffix)) is not None
                    and c not in (f"p9010{suffix}", f"p9050{suffix}", f"p5010{suffix}")]
        probs = np.array([p for _, p in pct_cols])
        order = np.argsort(probs)
        probs = probs[order]
        q = df[[pct_cols[i][0] for i in order]].iloc[0].to_numpy(dtype=float)
        mask = np.isfinite(q)
        probs, q = probs[mask], q[mask]
        from scipy.stats import norm
        sd = float(df[schema["d1_std"]].iloc[0])
        z = norm.ppf(probs)
        ax.plot(z, q, "o-", color=color, ms=4, label=f"{cc} {year} empirical")
        lims = ax.get_xlim()
        xx = np.linspace(min(z), max(z), 50)
        ax.plot(xx, xx * sd, "--", color="grey", lw=1,
                label=f"N(0,{sd:.2f}^2) reference")
        ax.set_xlabel("normal quantile")
        ax.set_ylabel("empirical quantile of 1y log change")
        ax.set_title(f"{cc} {year}")
        ax.legend(fontsize=7)

    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    pl_years = [2002, 2009, 2016]
    for ax, yr in zip(axes[0], pl_years):
        qq_panel(ax, "POL", yr, STYLE["POL"]["color"])
    qq_panel(axes[1, 0], "FRA", 2010, STYLE["FRA"]["color"])
    qq_panel(axes[1, 1], "USA", 2010, STYLE["USA"]["color"])
    axes[1, 2].axis("off")
    fig.suptitle("1-year log-change distribution vs normal (same variance):\n"
                 "a straight line on this plot = Gaussian; curvature at the ends = leptokurtosis")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(os.path.join(OUTDIR, "04_qq_shape.png"))
    plt.close(fig)

    return panels


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["A", "B", "all"], default="A")
    args = ap.parse_args()
    if args.stage in ("A", "all"):
        stage_a()
        print(f"Stage A done -> {OUTDIR}")
    if args.stage in ("B", "all"):
        print("Stage B not yet available from this invocation; run after the "
              "Poland estimate finishes (see scripts/pl_diagnostics_stageB.py).")
