"""
Sec 24.3: how much of each country's var_log_earns target is BETWEEN-AGE
variance (differences in mean log earnings across ages), which the model
cannot generate -- moments_from_panel computes the variance within each
age (panel column) and averages over ages.

For PL/FRA/GER/USA (Male, 25-55, same windows as the targets: 2001-2016 for
FRA/GER/USA, 2000-2020 for PL) computes the mean over years of
  std_resid^2   -- residualized log annual earnings (age dummies, per year)
  std_log_inc^2 -- raw log annual earnings (the current target concept)
and their difference (= between-age variance, to the extent both are on
the same sample -- checked via the observation counts).

Inputs: data/std resid.csv (std_resid, all four), data/Stats_20260720145608.csv
(FRA/GER/USA std_log_inc), data/PL std_log_inc.csv (PL std_loginc).
Writes output/resid_vs_loginc.json.
"""
from __future__ import annotations

import json

import pandas as pd

WINDOWS = {"POL": (2000, 2020), "FRA": (2001, 2016), "GER": (2001, 2016), "USA": (2001, 2016)}
RESID_CODE = {"POL": "POL", "FRA": "FRA", "GER": "DEU", "USA": "USA"}   # std resid.csv uses ISO DEU
TARGETS = {"POL": "targets/pol_grid_v2.json", "FRA": "targets/fra_grid.json",
           "GER": "targets/ger_grid.json", "USA": "targets/usa_grid.json"}


def sel(df, code, lo, hi):
    return df[(df.country == code) & (df.gender == "Male") & (df.age == "25-55")
              & df.year.between(lo, hi)].set_index("year").sort_index()


def main():
    resid = pd.read_csv("data/std resid.csv")
    stats = pd.read_csv("data/Stats_20260720145608.csv").rename(
        columns={"std_log_inc": "std_loginc", "nobs_log_inc": "nobs_loginc"})
    pl = pd.read_csv("data/PL std_log_inc.csv")
    pl_v1 = json.load(open("targets/pol_grid.json"))["var_log_earns"]

    out = {}
    for cc, (lo, hi) in WINDOWS.items():
        r = sel(resid, RESID_CODE[cc], lo, hi)
        l = sel(pl if cc == "POL" else stats, cc, lo, hi)
        assert list(r.index) == list(l.index) == list(range(lo, hi + 1)), cc
        v_res = (r.std_resid ** 2).mean()
        v_log = (l.std_loginc ** 2).mean()
        target = json.load(open(TARGETS[cc]))["var_log_earns"]
        nobs_ratio = (r.nobs_resid / l.nobs_loginc)
        out[cc] = {
            "window": f"{lo}-{hi}", "n_years": len(r),
            "var_resid": v_res, "var_loginc": v_log,
            "between_age": v_log - v_res,
            "between_age_share": (v_log - v_res) / v_log,
            "target_var_log_earns": target,
            "target_minus_var_loginc": target - v_log,
            "nobs_resid_over_nobs_loginc_min": float(nobs_ratio.min()),
            "nobs_resid_over_nobs_loginc_max": float(nobs_ratio.max()),
            "between_age_by_year": {int(y): float(l.std_loginc[y] ** 2 - r.std_resid[y] ** 2)
                                    for y in r.index},
        }
    out["POL"]["target_v1_std_perm"] = pl_v1

    print(f"{'cc':4s} {'window':10s} {'var_resid':>9s} {'var_loginc':>10s} {'between':>8s} "
          f"{'share':>6s} {'target':>7s} {'nobs r/l':>12s}")
    for cc, o in out.items():
        print(f"{cc:4s} {o['window']:10s} {o['var_resid']:9.4f} {o['var_loginc']:10.4f} "
              f"{o['between_age']:8.4f} {o['between_age_share']:6.1%} "
              f"{o['target_var_log_earns']:7.4f} "
              f"{o['nobs_resid_over_nobs_loginc_min']:.4f}-{o['nobs_resid_over_nobs_loginc_max']:.4f}")
    with open("output/resid_vs_loginc.json", "w") as f:
        json.dump(out, f, indent=2)


if __name__ == "__main__":
    main()
