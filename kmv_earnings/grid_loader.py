"""
Building the 8 target moments from GRID (Global Repository of Income
Dynamics, https://www.grid-database.org) downloads.

IMPORTANT / honest caveat: this file was written without live access to the
GRID download portal, so the exact CSV column names could not be verified.
The GRID country files follow a common master-code format (statistics by
year / gender / age group for: log earnings levels, 1-year and 5-year
residualized log earnings changes), but you WILL likely need to adjust
COLUMN_MAP below after downloading a file and inspecting its header.
Everything else in the pipeline is independent of this file: you can always
fall back to typing the 8 numbers into a targets JSON by hand.

What GRID gives vs what KMV needs:
  - Variance of log earnings: GRID reports sd / variance (and percentiles) of
    log earnings levels -> square the sd if needed.
  - Variance and kurtosis of 1y and 5y changes: reported directly
    (use the standard, moment-based kurtosis, NOT the Crow-Siddiqui robust
    kurtosis, to be comparable with KMV's 17.8 / 11.6).
  - "Frac 1yr change < X%": NOT reported directly by GRID. But GRID reports a
    fine set of percentiles of the change distribution; we invert the
    percentile function to get the CDF and read off
    P(|dlog y| < 0.1 / 0.2 / 0.5). See fractions_from_percentiles().

Comparability with KMV's US targets (GKOS SSA data): males, ~ages 25-55(60),
labor earnings above a minimum-earnings threshold, changes computed on
residualized log earnings (age effects removed). Use the GRID "male" x
prime-age cells and the residualized-change statistics, and average the
yearly statistics over the sample years.
"""

from __future__ import annotations

import re
import json
import numpy as np
import pandas as pd

# ---- adjust these after inspecting a real GRID csv header -------------------
COLUMN_MAP = {
    "year": "year",
    "gender": "gender",       # or "male" indicator / "sex"
    "age": "age",             # age group identifier
    "var_log": "var_logearn", # variance (or sd -> set SD_TO_VAR) of log earnings
    "var_d1": "var_d1_logearn",
    "var_d5": "var_d5_logearn",
    "kurt_d1": "kurt_d1_logearn",
    "kurt_d5": "kurt_d5_logearn",
}
SD_TO_VAR = False  # set True if the file reports sd instead of variance
# -----------------------------------------------------------------------------


def percentile_columns(df: pd.DataFrame, pattern: str = r"^p(\d{1,2}(?:_5)?)$"):
    """
    Find percentile columns like p1, p2_5, p5, ..., p99 and return
    (sorted probs in (0,1), matching column names).
    """
    probs, cols = [], []
    for c in df.columns:
        m = re.match(pattern, c)
        if m:
            probs.append(float(m.group(1).replace("_", ".")) / 100.0)
            cols.append(c)
    order = np.argsort(probs)
    return np.array(probs)[order], [cols[i] for i in order]


def fractions_from_percentiles(probs: np.ndarray, values: np.ndarray,
                               thresholds=(0.10, 0.20, 0.50)) -> dict:
    """
    Given percentiles (probs in (0,1), values = quantiles of dlog y),
    approximate P(|dlog y| < c) by interpolating the CDF:
        F(x) ~ interp of probs as a function of quantile values,
        P(|d| < c) = F(c) - F(-c).
    Accurate as long as the percentile grid is reasonably fine near +-c.
    """
    values = np.asarray(values, dtype=float)
    probs = np.asarray(probs, dtype=float)

    def F(x):
        return float(np.interp(x, values, probs, left=0.0, right=1.0))

    return {f"frac_d1_lt_{int(c*100)}": F(c) - F(-c) for c in thresholds}


def targets_from_grid_csvs(
    levels_csv: str,
    d1_csv: str,
    d5_csv: str,
    gender: str = "male",
    year_range: tuple[int, int] | None = None,
    column_map: dict | None = None,
) -> dict:
    """
    Build the 8-moment target dict from three GRID downloads:
      levels_csv - log earnings *levels* statistics
      d1_csv     - 1-year residualized log-change statistics (incl. percentiles)
      d5_csv     - 5-year residualized log-change statistics
    Yearly statistics are averaged over the selected years.
    """
    cm = {**COLUMN_MAP, **(column_map or {})}

    def prep(path):
        df = pd.read_csv(path)
        if cm["gender"] in df.columns:
            df = df[df[cm["gender"]].astype(str).str.lower().str.startswith(gender[0])]
        if year_range and cm["year"] in df.columns:
            df = df[df[cm["year"]].between(*year_range)]
        return df

    lv, d1, d5 = prep(levels_csv), prep(d1_csv), prep(d5_csv)

    var_log = lv[cm["var_log"]].mean()
    if SD_TO_VAR:
        var_log = (lv[cm["var_log"]] ** 2).mean()

    targets = {
        "var_log_earns": float(var_log),
        "var_d1": float(d1[cm["var_d1"]].mean()),
        "var_d5": float(d5[cm["var_d5"]].mean()),
        "kurt_d1": float(d1[cm["kurt_d1"]].mean()),
        "kurt_d5": float(d5[cm["kurt_d5"]].mean()),
    }

    # fractions of small changes from the percentile grid of the 1y change file
    probs, cols = percentile_columns(d1)
    if len(cols) >= 5:
        q = d1[cols].mean(axis=0).to_numpy()
        targets.update(fractions_from_percentiles(probs, q))
    else:
        raise ValueError(
            "No percentile columns found in the 1-year change file; "
            "adjust percentile_columns() pattern or enter the fractions manually."
        )
    return targets


def save_targets(targets: dict, path: str) -> None:
    with open(path, "w") as f:
        json.dump(targets, f, indent=2)


def load_targets(path: str) -> dict:
    with open(path) as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Loader for GRID "Stats_*.csv"-style exports: one wide file, rows = country x
# year, columns giving level dispersion, 1y/5y residualized log-change
# dispersion, and a percentile grid of the 1y change. GRID's own column naming
# has changed between downloads (verified on two real exports, July 2026 and
# October 2026) -- the schema below is auto-detected per file, not hardcoded,
# so a third naming convention just needs a new entry in _SCHEMAS.
#
#   schema "res_1yr_log_chg" (FRA/GER/USA download, July 2026):
#     nobs_log_inc, std_log_inc,
#     nobs_res_1yr_log_chg, std_res_1yr_log_chg, kurt_res_1yr_log_chg,
#     p{q}_res_1yr_log_chg for q in 1, 2_5, 5, 10, 12_5, 25, ..., 99_9, 99_99,
#     nobs_res_5yr_log_chg, std_res_5yr_log_chg, kurt_res_5yr_log_chg
#
#   schema "resid_dNyr" (Poland download, October 2026):
#     nobs_perm, std_perm,
#     nobs_resid_d1yr, std_resid_d1yr, kurt_resid_d1yr,
#     p{q}_resid_d1yr for the same q grid,
#     nobs_resid_d5yr, std_resid_d5yr, kurt_resid_d5yr
# ---------------------------------------------------------------------------

_SCHEMAS = [
    {
        "name": "res_1yr_log_chg",
        "level_std": "std_log_inc",
        "d1_suffix": "_res_1yr_log_chg",
        "d1_nobs": "nobs_res_1yr_log_chg",
        "d1_std": "std_res_1yr_log_chg",
        "d1_kurt": "kurt_res_1yr_log_chg",
        "d5_nobs": "nobs_res_5yr_log_chg",
        "d5_std": "std_res_5yr_log_chg",
        "d5_kurt": "kurt_res_5yr_log_chg",
    },
    {
        "name": "resid_dNyr",
        "level_std": "std_perm",
        "d1_suffix": "_resid_d1yr",
        "d1_nobs": "nobs_resid_d1yr",
        "d1_std": "std_resid_d1yr",
        "d1_kurt": "kurt_resid_d1yr",
        "d5_nobs": "nobs_resid_d5yr",
        "d5_std": "std_resid_d5yr",
        "d5_kurt": "kurt_resid_d5yr",
    },
]


def _detect_schema(df: pd.DataFrame) -> dict:
    for schema in _SCHEMAS:
        if schema["level_std"] in df.columns and schema["d1_std"] in df.columns:
            return schema
    raise ValueError(
        "Unrecognised GRID stats column schema -- none of the known column "
        f"sets match {sorted(df.columns)}. Add a new entry to _SCHEMAS in "
        "grid_loader.py after inspecting the file's header."
    )


def _pct_prob_from_name(col: str, suffix: str) -> float | None:
    """'p12_5_res_1yr_log_chg' -> 0.125, 'p99_99_...' -> 0.9999."""
    if not (col.startswith("p") and col.endswith(suffix)):
        return None
    body = col[1:-len(suffix)].rstrip("_")
    if not body or not body[0].isdigit():
        return None
    return float(body.replace("_", ".")) / 100.0


def targets_from_grid_stats_csv(
    path: str,
    country: str,
    year_range: tuple[int, int] | None = None,
    weight_by_nobs: bool = False,
    gender: str | None = "Male",
) -> dict:
    """
    Build the 8-moment target dict for one country from a GRID Stats_*.csv.
    Yearly statistics are averaged over the selected years (simple mean by
    default; weight_by_nobs=True weights by cell observation counts).
    Fractions of small 1y changes are interpolated from the percentile grid
    of the residualized 1-year change distribution.

    `gender`: filter the `gender` column to this value before anything else
    (default "Male", matching the KMV/GKOS sample design, §10). Pass None to
    skip filtering -- but then an assertion below still requires that exactly
    one gender value survives, because averaging statistics across genders
    silently mixes in participation-margin effects that have nothing to do
    with the income *process* (§10). Raise, don't average, if that is not
    already the case.
    """
    df = pd.read_csv(path)
    df = df[df["country"] == country]
    if gender is not None and "gender" in df.columns:
        df = df[df["gender"].astype(str).str.lower() == gender.lower()]
    if year_range:
        df = df[df["year"].between(*year_range)]
    if df.empty:
        raise ValueError(f"No rows for country={country} in {path}")

    if "gender" in df.columns:
        genders = df["gender"].dropna().unique()
        assert len(genders) == 1, (
            f"targets_from_grid_stats_csv would silently average over "
            f"genders {sorted(map(str, genders))} for country={country!r} -- "
            "pass gender='Male' (or whichever single value the file uses) to "
            "filter first."
        )

    schema = _detect_schema(df)

    def avg(col, wcol=None):
        if weight_by_nobs and wcol is not None:
            w = df[wcol]
            return float((df[col] * w).sum() / w.sum())
        return float(df[col].mean())

    df = df.assign(_v_log=df[schema["level_std"]] ** 2,
                   _v_d1=df[schema["d1_std"]] ** 2,
                   _v_d5=df[schema["d5_std"]] ** 2)
    targets = {
        "var_log_earns": avg("_v_log"),
        "var_d1": avg("_v_d1", schema["d1_nobs"]),
        "var_d5": avg("_v_d5", schema["d5_nobs"]),
        "kurt_d1": avg(schema["d1_kurt"], schema["d1_nobs"]),
        "kurt_d5": avg(schema["d5_kurt"], schema["d5_nobs"]),
    }

    suffix = schema["d1_suffix"]
    probs, cols = [], []
    for c in df.columns:
        p = _pct_prob_from_name(c, suffix)
        if p is not None and c not in ("p9010" + suffix, "p9050" + suffix,
                                       "p5010" + suffix):
            probs.append(p)
            cols.append(c)
    order = np.argsort(probs)
    probs = np.array(probs)[order]
    q = df[[cols[i] for i in order]].mean(axis=0).to_numpy()
    targets.update(fractions_from_percentiles(probs, q))
    return targets
