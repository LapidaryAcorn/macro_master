"""
Sec 26 (A.1, A.3, A.4): Poland v3 targets and linear trends in the moments.

v3 = v2 with the level moment residualized: var_log_earns = mean over
2000-2020 of GRID std_resid^2 (data/std resid.csv, POL, Male, 25-55). The other
seven moments are the v2 ones (already residualized changes), copied from
targets/pol_grid_v2.json.

Year panel (output/pol_v3/year_panel_pol_v3.csv) = output/pl_diagnostics/
year_panel_pol_v2.csv with var_log replaced by std_resid^2. GRID indexes a
change by its BASE year: d1 in year t is t -> t+1 (series 2000-2019), d5 in
year t is t -> t+5 (series 2000-2015); the level runs 2000-2020.

Trends: OLS of each moment's annual series on a constant and the year, using
each series' OWN available years (no padding, no common sub-window). Slope SE
is Newey-West (Bartlett kernel, lag 2 = floor(4 (T/100)^(2/9)) for T = 16..21,
no small-sample correction). Trend at the window midpoint (mean of the series'
own years) equals the series mean exactly (OLS with a constant) -- checked.

v3_end (A.4): trend value at base year 2019 (last year with 1-year changes)
for every moment; for var_d5 / kurt_d5 that is an extrapolation 4 years past
the last data year (2015). Robustness test only, not the baseline.

Outputs: targets/pol_grid_v3.json, targets/pol_grid_v3_end.json,
output/pol_v3/year_panel_pol_v3.csv, output/pol_v3/trends/{trends.csv,
trends.json, trends.png, trend_residuals.csv}.
"""

from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from kmv_earnings.simulate import MOMENT_ORDER, MOMENT_LABELS

WINDOW = (2000, 2020)
END_YEAR = 2019            # last base year with a 1-year change
NW_LAG = 2
OUT = "output/pol_v3"
TR = os.path.join(OUT, "trends")

COLS = {  # target key -> year-panel column
    "var_log_earns": "var_log", "var_d1": "var_d1", "var_d5": "var_d5",
    "kurt_d1": "kurt_d1", "kurt_d5": "kurt_d5",
    "frac_d1_lt_10": "frac_d1_lt_10", "frac_d1_lt_20": "frac_d1_lt_20",
    "frac_d1_lt_50": "frac_d1_lt_50",
}


def build_panel() -> pd.DataFrame:
    v2 = pd.read_csv("output/pl_diagnostics/year_panel_pol_v2.csv")
    sr = pd.read_csv("data/std resid.csv")
    sr = sr[(sr.country == "POL") & (sr.gender == "Male") & (sr.age == "25-55")]
    sr = sr.set_index("year")
    p = v2.copy()
    p["var_log_v2"] = p["var_log"]
    p["var_log"] = p["year"].map(sr["std_resid"] ** 2)
    p["nobs_resid"] = p["year"].map(sr["nobs_resid"])
    assert p["var_log"].notna().all()
    assert (p["nobs_resid"] == p["nobs_level"]).all()  # same sample [M], Sec 24.3
    return p


def newey_west_ols(t: np.ndarray, y: np.ndarray, lag: int):
    X = np.column_stack([np.ones_like(t), t])
    XtX_inv = np.linalg.inv(X.T @ X)
    b = XtX_inv @ X.T @ y
    u = y - X @ b
    Xu = X * u[:, None]
    S = Xu.T @ Xu
    for l in range(1, lag + 1):
        w = 1.0 - l / (lag + 1.0)
        G = Xu[l:].T @ Xu[:-l]
        S += w * (G + G.T)
    V = XtX_inv @ S @ XtX_inv
    se_ols = np.sqrt(np.diag(XtX_inv) * (u @ u) / (len(y) - 2))
    return b, np.sqrt(np.diag(V)), se_ols, u


def main():
    os.makedirs(TR, exist_ok=True)
    panel = build_panel()
    panel = panel[(panel.year >= WINDOW[0]) & (panel.year <= WINDOW[1])]
    panel.to_csv(os.path.join(OUT, "year_panel_pol_v3.csv"), index=False)

    v2 = json.load(open("targets/pol_grid_v2.json"))
    v3 = dict(v2)
    v3["var_log_earns"] = float(panel["var_log"].mean())
    # variances/kurtoses are reproduced from the panel exactly. The three
    # fractions are not: the targets interpolate the CDF from the across-year
    # MEAN percentiles (grid_loader), the panel interpolates year by year, so
    # the panel mean differs by <= 0.006 (Sec 26.1). Targets stay as in v2.
    for k in MOMENT_ORDER[1:]:
        tol = 1e-12 if not k.startswith("frac") else 0.01
        assert abs(panel[COLS[k]].mean() - v2[k]) < tol, k
    with open("targets/pol_grid_v3.json", "w") as f:
        json.dump(v3, f, indent=2)

    rows, resid = [], pd.DataFrame({"year": panel.year.values})
    v3_end, end_meta = {}, {}
    for k in MOMENT_ORDER:
        s = panel[["year", COLS[k]]].dropna()
        t, y = s.year.to_numpy(float), s[COLS[k]].to_numpy(float)
        tmid = t.mean()
        b, se_nw, se_ols, u = newey_west_ols(t - tmid, y, NW_LAG)
        trend = lambda yr: float(b[0] + b[1] * (yr - tmid))
        assert abs(trend(tmid) - y.mean()) < 1e-12
        last = int(t.max())
        rows.append({
            "moment": k, "label": MOMENT_LABELS[k],
            "years": f"{int(t.min())}-{last}", "T": len(t),
            "mean": y.mean(), "slope_per_yr": b[1], "se_nw": se_nw[1],
            "t_nw": b[1] / se_nw[1], "se_ols": se_ols[1],
            "trend_2000": trend(2000), "mid_year": tmid, "trend_mid": trend(tmid),
            "trend_end": trend(last), "end_year": last,
            "trend_2019": trend(END_YEAR),
            "extrap_years_2019": max(0, END_YEAR - last),
            "slope_x_window_rel": b[1] * (t.max() - t.min()) / y.mean(),
        })
        resid[k] = pd.Series(dict(zip(t.astype(int), u))).reindex(resid.year).values
        # v3 target + trend shift, so v3_end - v3 is the trend alone (for the
        # fractions this keeps the target's mean-percentile construction)
        v3_end[k] = v3[k] + float(b[1]) * (END_YEAR - tmid)
        end_meta[k] = {"data_years": f"{int(t.min())}-{last}",
                       "extrapolated": END_YEAR > last,
                       "extrapolation_years": max(0, END_YEAR - last)}

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(TR, "trends.csv"), index=False)
    resid.to_csv(os.path.join(TR, "trend_residuals.csv"), index=False)
    with open(os.path.join(TR, "trends.json"), "w") as f:
        json.dump({"nw_lag": NW_LAG, "rows": rows}, f, indent=2, default=float)
    with open("targets/pol_grid_v3_end.json", "w") as f:
        json.dump(v3_end, f, indent=2)
    with open(os.path.join(TR, "v3_end_meta.json"), "w") as f:
        json.dump({"end_base_year": END_YEAR, "moments": end_meta,
                   "note": "trend value at GRID base year 2019; d5 series end "
                           "in 2015 (2015->2020), so their 2019 value is a "
                           "4-year linear extrapolation"}, f, indent=2)

    pd.set_option("display.width", 200)
    print(df[["moment", "years", "mean", "slope_per_yr", "se_nw", "t_nw",
              "trend_2000", "trend_mid", "trend_end", "trend_2019",
              "extrap_years_2019"]].to_string(index=False, float_format="%.4g"))
    print("\nv2 / v3 / v3_end targets:")
    for k in MOMENT_ORDER:
        print(f"  {k:15s} {v2[k]:9.4f} {v3[k]:9.4f} {v3_end[k]:9.4f}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker
    fig, axes = plt.subplots(2, 4, figsize=(16, 7))
    for ax, r in zip(axes.ravel(), rows):
        k = r["moment"]
        s = panel[["year", COLS[k]]].dropna()
        ax.plot(s.year, s[COLS[k]], "o-", ms=3, lw=1, color="C0", label="GRID, annual")
        yr = np.arange(2000, END_YEAR + 1)
        tr = r["trend_mid"] + r["slope_per_yr"] * (yr - r["mid_year"])
        inside = yr <= r["end_year"]
        ax.plot(yr[inside], tr[inside], "-", color="C3", label="linear trend (OLS)")
        if (~inside).any():
            ax.plot(yr[r["end_year"] - 2000:], tr[r["end_year"] - 2000:], ":",
                    color="C3", label="extrapolation")
        ax.axhline(r["mean"], color="0.6", lw=0.8, ls="--", label="window mean (v3)")
        ax.plot([END_YEAR], [r["trend_2019"]], "s", color="C3", ms=5)
        ax.set_title(f"{MOMENT_LABELS[k]}\nslope {r['slope_per_yr']:+.3g}/yr "
                     f"(NW t = {r['t_nw']:.1f})", fontsize=9)
        ax.tick_params(labelsize=8)
        ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(5))
    axes[0, 0].legend(fontsize=7)
    fig.suptitle("Poland (GRID, men 25-55): annual moments and linear trends; "
                 "square = v3_end target (trend at base year 2019)", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(TR, "trends.png"), dpi=130)
    print(f"\nSaved: targets/pol_grid_v3.json, targets/pol_grid_v3_end.json, {TR}/")


if __name__ == "__main__":
    main()
