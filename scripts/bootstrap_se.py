"""
Task 3 (DECISIONS.md, 2026-10-11): bootstrap standard errors for the KMV
earnings-process estimates, since GRID does not publish SE on its reported
statistics but does publish them year by year.

Design
------
For each country: resample years WITH REPLACEMENT from that country's own
estimation window (same number of draws as years in the window), recompute
the 8 target moments as the simple across-year mean over the resampled years
(verified to reproduce the original targets/<cc>_grid.json bit-for-bit when
applied to the *unresampled* window -- see scripts/pl_diagnostics.py's
year_panel_*.csv, which this script reuses), then re-estimate the KMV process
with a WARM START at the country's own point estimate and a SHRUNK optimizer
budget (reduced simulated-panel size, DE iterations/population, and NM polish
iterations) -- re-running the full serious budget >=200 times per country
would cost tens of hours (see the time-estimate given to the user before this
script was run).

Country restrictions (fixed/derive) are kept IDENTICAL to the canonical
estimation (scripts/close_step1.py): GER and USA still pin beta2 via
ErgodicVarConstraint in every bootstrap rep, so their reported SE/CI only
cover the five free parameters -- beta2 there is calibrated, not estimated,
and bootstrapping a calibrated quantity would misrepresent it as sampling
uncertainty.

CRN note: model_moments() uses a FIXED internal simulation seed by default
(DECISIONS Sec 6). This script does not override it, so every bootstrap rep's
re-estimation shares the same simulated-worker draws; the only source of
variation across reps is which GRID years were resampled into the target
moments. This is intentional -- it isolates year-to-year sampling variability
in the DATA from simulation noise, which is a separate (and, by Sec 6,
already-controlled) concern.

Budget-validation sub-step (required by the task): before committing to the
reduced budget, N_CHECK resamples per country are also re-estimated with a
larger REFERENCE budget, and the two are compared parameter-by-parameter. If
the reduced budget shows a systematic (not just noisy) shift relative to the
reference, STOP and widen the reduced budget before running the full sweep.

Honest caveat (stated in the task): this measures year-to-year *sampling*
variability in the reported GRID statistics, not full sampling uncertainty
(which would require microdata and is not available here).
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from kmv_earnings.simulate import PARAM_ORDER, MOMENT_ORDER
from kmv_earnings.estimate import estimate, ErgodicVarConstraint, load_params

KMV_ERGODIC_VAR2 = 0.007 * 1.53 ** 2 / (2.0 * 0.009)
_PIN_BETA2 = {"beta2": ErgodicVarConstraint("2", target_var=KMV_ERGODIC_VAR2)}

COUNTRIES = {
    "USA": dict(targets="targets/usa_grid.json", params="output/usa/params.json",
                panel="output/pl_diagnostics/year_panel_usa.csv", window=(2001, 2016),
                derive=dict(_PIN_BETA2)),
    "FRA": dict(targets="targets/fra_grid.json", params="output/fra/params.json",
                panel="output/pl_diagnostics/year_panel_fra.csv", window=(2001, 2016),
                derive={}),
    "GER": dict(targets="targets/ger_grid.json", params="output/ger/params.json",
                panel="output/pl_diagnostics/year_panel_ger.csv", window=(2001, 2016),
                derive=dict(_PIN_BETA2)),
    "POL": dict(targets="targets/pol_grid.json", params="output/pol/params.json",
                panel="output/pl_diagnostics/year_panel_pol.csv", window=(2000, 2020),
                derive={}),
}

MOMENT_COLS = {
    "var_log_earns": "var_log", "var_d1": "var_d1", "var_d5": "var_d5",
    "kurt_d1": "kurt_d1", "kurt_d5": "kurt_d5",
    "frac_d1_lt_10": "frac_d1_lt_10", "frac_d1_lt_20": "frac_d1_lt_20",
    "frac_d1_lt_50": "frac_d1_lt_50",
}

REDUCED = dict(n_workers_sim=5_000, de_maxiter=10, de_popsize=5, nm_maxiter=50)
REFERENCE = dict(n_workers_sim=15_000, de_maxiter=20, de_popsize=6, nm_maxiter=100)

OUTDIR = "output/bootstrap"

# Sec 24: --version v2 swaps Poland to the v2 targets / point estimate and a
# year panel whose var_log is GRID std_loginc^2 (output/pl_diagnostics/
# year_panel_pol_v2.csv; v1 panel untouched), writes to output/bootstrap_v2/.
V2_OVERRIDES = {"POL": dict(targets="targets/pol_grid_v2.json",
                            params="output/pol_v2/params.json",
                            panel="output/pl_diagnostics/year_panel_pol_v2.csv")}


def _apply_version(version: str):
    # Pool workers on Windows re-import this module (spawn), so the version is
    # passed through the environment and applied at import time in every
    # process, not only in main().
    global OUTDIR
    if version == "v2":
        for cc, ov in V2_OVERRIDES.items():
            COUNTRIES[cc].update(ov)
        OUTDIR = "output/bootstrap_v2"


_apply_version(os.environ.get("KMV_BOOTSTRAP_VERSION", "v1"))


def resample_targets(panel: pd.DataFrame, window: tuple, rng: np.random.Generator) -> dict:
    yrs = panel[(panel.year >= window[0]) & (panel.year <= window[1])]["year"].to_numpy()
    draw = rng.choice(yrs, size=len(yrs), replace=True)
    sub = panel.set_index("year").loc[draw]
    return {k: float(np.nanmean(sub[c].to_numpy())) for k, c in MOMENT_COLS.items()}


def _one_rep(args):
    cc, rep_seed, budget_name = args
    cfg = COUNTRIES[cc]
    panel = pd.read_csv(cfg["panel"])
    rng = np.random.default_rng(rep_seed)
    targets = resample_targets(panel, cfg["window"], rng)
    budget = REDUCED if budget_name == "reduced" else REFERENCE
    x0 = load_params(cfg["params"])
    t0 = time.time()
    params, res = estimate(
        targets,
        sim_kwargs=dict(n_workers=budget["n_workers_sim"], n_years_keep=36),
        maxiter=budget["de_maxiter"], popsize=budget["de_popsize"],
        workers=1, seed=rep_seed, x0=x0, disp=False, polish_nm=True,
        nm_maxiter=budget["nm_maxiter"], derive=cfg["derive"],
    )
    return {"country": cc, "rep_seed": rep_seed, "budget": budget_name,
            "params": params, "objective": float(res.fun),
            "wall_s": time.time() - t0, "targets": targets}


def run_validation(n_check: int, pool: Pool) -> list[dict]:
    jobs = []
    for cc_idx, cc in enumerate(COUNTRIES):
        for i in range(n_check):
            seed = 9_000_000 + cc_idx * 100 + i
            jobs.append((cc, seed, "reduced"))
            jobs.append((cc, seed, "reference"))
    print(f"[validation] {len(jobs)} calls ({n_check} resamples x {len(COUNTRIES)} "
          "countries x 2 budgets) ...", flush=True)
    results = pool.map(_one_rep, jobs)
    with open(os.path.join(OUTDIR, "validation.json"), "w") as f:
        json.dump(results, f, indent=2, default=float)
    return results


def report_validation(results: list[dict]):
    by_key = {}
    for r in results:
        by_key.setdefault((r["country"], r["rep_seed"]), {})[r["budget"]] = r
    print("\n[validation] reduced vs reference budget, same resampled targets:")
    print(f"  {'country':5s} {'seed':>10s} {'param':8s} {'reduced':>10s} "
          f"{'reference':>10s} {'rel.diff':>9s}")
    by_cc = {}
    for (cc, seed), pair in by_key.items():
        if "reduced" not in pair or "reference" not in pair:
            continue
        for p in PARAM_ORDER:
            a, b = pair["reduced"]["params"][p], pair["reference"]["params"][p]
            rd = (a - b) / max(abs(a), abs(b), 1e-12)
            by_cc.setdefault(cc, {}).setdefault(p, []).append(rd)
            print(f"  {cc:5s} {seed:>10d} {p:8s} {a:10.5g} {b:10.5g} {rd:>+8.1%}")
    print("\n[validation] summary PER COUNTRY per parameter (mean rel.diff +- sd; "
          "checked per country, not pooled, so one bad country can't hide behind "
          "the others' good behaviour):")
    verdict_ok = True
    for cc, pdict in by_cc.items():
        for p in PARAM_ORDER:
            d = np.array(pdict.get(p, []))
            if len(d) == 0:
                continue
            mean, sd = d.mean(), d.std()
            flag = ""
            if abs(mean) > 0.15 and (len(d) < 3 or abs(mean) > sd):
                flag = "  <-- possible systematic shift"
                verdict_ok = False
            print(f"  {cc:5s} {p:8s} mean={mean:+.1%}  sd={sd:.1%}  n={len(d)}{flag}")
    print(f"\n[validation] verdict: {'OK, reduced budget is not systematically biased' if verdict_ok else 'WARNING -- widen REDUCED budget (or check for a stale warm start) before trusting the full sweep'}")
    return verdict_ok


def run_bootstrap(n_boot: int, pool: Pool):
    jobs = []
    for cc in COUNTRIES:
        for i in range(n_boot):
            jobs.append((cc, 1_000_000 + i, "reduced"))
    print(f"\n[bootstrap] {len(jobs)} reps ({n_boot} per country x {len(COUNTRIES)} "
          "countries), reduced budget, warm-started ...", flush=True)
    t0 = time.time()
    results = []
    chunk = 200
    for i in range(0, len(jobs), chunk):
        batch = jobs[i:i + chunk]
        results.extend(pool.map(_one_rep, batch))
        elapsed = time.time() - t0
        print(f"  {len(results)}/{len(jobs)} done, {elapsed:.0f}s elapsed, "
              f"~{elapsed/max(len(results),1)*(len(jobs)-len(results)):.0f}s remaining",
              flush=True)
    with open(os.path.join(OUTDIR, "bootstrap_reps.json"), "w") as f:
        json.dump(results, f, indent=2, default=float)
    summarize(results)


def summarize(results: list[dict]):
    rows = []
    for cc, cfg in COUNTRIES.items():
        reps = [r for r in results if r["country"] == cc]
        if not reps:
            continue
        locked = set(cfg["derive"])
        row = {"country": cc, "n_reps": len(reps)}
        for p in PARAM_ORDER:
            vals = np.array([r["params"][p] for r in reps])
            row[f"{p}_mean"] = float(vals.mean())
            row[f"{p}_se"] = float(vals.std(ddof=1)) if p not in locked else None
            row[f"{p}_p2.5"] = float(np.percentile(vals, 2.5)) if p not in locked else None
            row[f"{p}_p97.5"] = float(np.percentile(vals, 97.5)) if p not in locked else None
            row[f"{p}_locked"] = p in locked
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUTDIR, "bootstrap_summary.csv"), index=False)
    print("\n[bootstrap] SE and 95% percentile CIs (locked params show calibrated value, no SE):")
    for row in rows:
        print(f"\n  {row['country']} (n={row['n_reps']}):")
        for p in PARAM_ORDER:
            if row[f"{p}_locked"]:
                print(f"    {p:8s} = {row[f'{p}_mean']:.5g}  [calibrated/derived, not bootstrapped]")
            else:
                print(f"    {p:8s} mean={row[f'{p}_mean']:.5g}  SE={row[f'{p}_se']:.3g}  "
                      f"95% CI=[{row[f'{p}_p2.5']:.5g}, {row[f'{p}_p97.5']:.5g}]")
    print(f"\nSaved: {OUTDIR}/bootstrap_summary.csv, bootstrap_reps.json")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-boot", type=int, default=200)
    ap.add_argument("--n-check", type=int, default=3)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument("--skip-validation", action="store_true")
    ap.add_argument("--version", choices=["v1", "v2"], default="v1")
    ap.add_argument("--countries", nargs="+", default=None)
    args = ap.parse_args()

    os.environ["KMV_BOOTSTRAP_VERSION"] = args.version
    _apply_version(args.version)
    if args.countries:
        for cc in list(COUNTRIES):
            if cc not in args.countries:
                del COUNTRIES[cc]
    os.makedirs(OUTDIR, exist_ok=True)
    with Pool(args.workers) as pool:
        if not args.skip_validation:
            ok = report_validation(run_validation(args.n_check, pool))
            if not ok:
                print("\nABORTING full bootstrap sweep -- see warning above. "
                      "Widen REDUCED in this file and re-run.")
                return
        run_bootstrap(args.n_boot, pool)


if __name__ == "__main__":
    main()
