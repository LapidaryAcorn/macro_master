"""
Entry point. Usage:

  # sanity check: simulate at KMV's published parameters, compare with Table 2
  python -m kmv_earnings.run validate

  # full SMM estimation against a targets file, then produce the table
  python -m kmv_earnings.run estimate --targets targets/us_kmv.json --out output/us

  # rebuild the table + chain export from saved parameters (no re-estimation)
  python -m kmv_earnings.run table --targets targets/us_kmv.json \
      --params output/us/params.json --out output/us

Flags for estimate: --quick (small simulation + few iterations, pipeline test
only), --workers N (parallel differential evolution).
"""

from __future__ import annotations

import argparse
import os
import pandas as pd

from .simulate import (MOMENT_ORDER, MOMENT_LABELS, KMV_TABLE3_PARAMS,
                       KMV_US_DATA_TARGETS, model_moments)
from .estimate import estimate, save_params, load_params
from .tiktak import tiktak
from .discretize import discretize_process, discretized_moments, export_chain
from .grid_loader import load_targets, save_targets


def build_table(targets: dict, params: dict, sim_kwargs: dict,
                disc_kwargs: dict | None = None) -> pd.DataFrame:
    m_est = model_moments(params, **sim_kwargs)
    disc = discretize_process(params, **(disc_kwargs or {}))
    disc_keys = {"n_workers", "n_years_keep", "lifecycle"}
    m_disc = discretized_moments(params, disc=disc,
                                 **{k: v for k, v in sim_kwargs.items() if k in disc_keys})
    rows = []
    for key in MOMENT_ORDER:
        rows.append({
            "Moment": MOMENT_LABELS[key],
            "Data": targets[key],
            "Model Estimated": round(m_est[key], 3),
            "Model Discretized": round(m_disc[key], 3),
        })
    return pd.DataFrame(rows), disc


def cmd_validate(args):
    sim_kwargs = dict(n_workers=args.n_workers, n_years_keep=36)
    print("Simulating at KMV (2018) Table 3 parameters...")
    df, disc = build_table(KMV_US_DATA_TARGETS, KMV_TABLE3_PARAMS, sim_kwargs)
    print(df.to_string(index=False))
    print("\n('Data' = KMV Table 2 US targets. 'Model Estimated' should be close to")
    print(" KMV's model column: 0.70 0.23 0.46 16.5 12.1 0.56 0.67 0.85.)")


def cmd_estimate(args):
    targets = load_targets(args.targets)
    if args.quick:
        sim_kwargs = dict(n_workers=10_000, n_years_keep=36, steps_per_quarter=3)
        maxiter, popsize = args.maxiter or 8, 6
    else:
        sim_kwargs = dict(n_workers=args.n_workers, n_years_keep=36)
        maxiter, popsize = args.maxiter or 100, 12

    if args.method == "tiktak":
        params, _ = tiktak(targets, sim_kwargs=sim_kwargs,
                           n_sobol=args.n_sobol, n_local=args.n_local)
    else:
        x0 = KMV_TABLE3_PARAMS if args.warm_start else None
        params, res = estimate(targets, sim_kwargs=sim_kwargs, maxiter=maxiter,
                               popsize=popsize, workers=args.workers, x0=x0)
    os.makedirs(args.out, exist_ok=True)
    save_params(params, os.path.join(args.out, "params.json"))
    print("\nEstimated parameters (quarterly rates):")
    for k, v in params.items():
        print(f"  {k:8s} = {v:.4f}")
    _make_outputs(targets, params, sim_kwargs, args.out)


def cmd_table(args):
    targets = load_targets(args.targets)
    params = load_params(args.params)
    sim_kwargs = dict(n_workers=args.n_workers, n_years_keep=36)
    os.makedirs(args.out, exist_ok=True)
    _make_outputs(targets, params, sim_kwargs, args.out)


def _make_outputs(targets, params, sim_kwargs, outdir):
    df, disc = build_table(targets, params, sim_kwargs)
    print("\n" + df.to_string(index=False))
    df.to_csv(os.path.join(outdir, "table2.csv"), index=False)
    with open(os.path.join(outdir, "table2.tex"), "w") as f:
        f.write(df.to_latex(index=False, float_format="%.2f",
                            caption="Earnings Process Estimation Fit",
                            label="tab:earnings_fit"))
    export_chain(disc, outdir)
    print(f"\nSaved: table2.csv, table2.tex, income_process_*.txt in {outdir}/")


def main():
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)

    pv = sub.add_parser("validate")
    pv.add_argument("--n-workers", type=int, default=50_000)
    pv.set_defaults(func=cmd_validate)

    pe = sub.add_parser("estimate")
    pe.add_argument("--targets", required=True)
    pe.add_argument("--out", required=True)
    pe.add_argument("--quick", action="store_true")
    pe.add_argument("--method", choices=["de", "tiktak"], default="de")
    pe.add_argument("--n-sobol", type=int, default=256)
    pe.add_argument("--n-local", type=int, default=10)
    pe.add_argument("--warm-start", action="store_true",
                    help="seed the population with the KMV US estimates")
    pe.add_argument("--maxiter", type=int, default=None)
    pe.add_argument("--workers", type=int, default=1)
    pe.add_argument("--n-workers", type=int, default=50_000)
    pe.set_defaults(func=cmd_estimate)

    pt = sub.add_parser("table")
    pt.add_argument("--targets", required=True)
    pt.add_argument("--params", required=True)
    pt.add_argument("--out", required=True)
    pt.add_argument("--n-workers", type=int, default=50_000)
    pt.set_defaults(func=cmd_table)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
