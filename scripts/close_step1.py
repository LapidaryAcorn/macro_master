"""
Close-out driver for step 1.

For each country (USA, FRA, GER):
  1. Full-budget SMM estimation with BOTH global optimisers (differential
     evolution and TikTak), serious simulation settings, warm-started at the
     KMV Table 3 parameters.
  2. Out-of-sample scoring of both estimates on a fresh-seed, larger panel.
  3. Agreement report: persistent-component parameters (lambda2/beta2/sigma2),
     bound-proximity of every parameter, and fit quality (8-moment table).
  4. The better estimate (lower OOS objective) becomes the canonical one:
     regenerate table2.{csv,tex} and the exported Markov chain
     income_process_*.txt, plus discretization diagnostics.

Targets are the common-window (2001-2016) GRID targets in targets/.

Usage (from the repo root):
    python scripts/close_step1.py                 # all three countries, serious
    python scripts/close_step1.py --countries USA
    python scripts/close_step1.py --quick         # pipeline smoke test
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from kmv_earnings.simulate import (MOMENT_ORDER, MOMENT_LABELS, PARAM_ORDER,
                                   KMV_TABLE3_PARAMS, model_moments)
from kmv_earnings.estimate import (estimate, DEFAULT_BOUNDS, save_params,
                                   ErgodicVarConstraint)
from kmv_earnings.tiktak import tiktak
from kmv_earnings.discretize import (discretize_process, discretized_moments,
                                     diagnose, export_chain,
                                     theoretical_stationary_var)
from kmv_earnings.run import build_table

COUNTRIES = {
    "USA": "targets/usa_grid.json",
    "FRA": "targets/fra_grid.json",
    "GER": "targets/ger_grid.json",
}
PERSISTENT = ("lambda2", "beta2", "sigma2")

# KMV (2018) Table 3 ergodic variance of the persistent component,
# lambda2 * sigma2^2 / (2 beta2) at (0.007, 0.009, 1.53).
KMV_ERGODIC_VAR2 = 0.007 * 1.53 ** 2 / (2.0 * 0.009)

# Per-country parameter restrictions (see close-out decisions):
#   GER's persistent-component mean reversion (beta2) is not identified by the
#   16-year German panel; pin it so the ergodic variance of that component
#   equals KMV's, and estimate the other five parameters.
COUNTRY_DERIVE = {
    "GER": {"beta2": ErgodicVarConstraint("2", target_var=KMV_ERGODIC_VAR2)},
}
COUNTRY_FIXED: dict[str, dict] = {}


def rel_objective(m: dict, targets: dict) -> float:
    t = np.array([targets[k] for k in MOMENT_ORDER])
    mv = np.array([m[k] for k in MOMENT_ORDER])
    return float(np.sum(((mv - t) / t) ** 2))


def oos_score(params: dict, targets: dict, sim_kwargs: dict, seed: int):
    """Score a parameter vector on a fresh-seed, larger panel (out of sample
    relative to the common-random-numbers seed used during optimisation)."""
    sk = dict(sim_kwargs)
    sk["n_workers"] = max(sk.get("n_workers", 50_000), 100_000)
    sk["seed"] = seed
    m = model_moments(params, **sk)
    return rel_objective(m, targets), m


def bound_table(params: dict, locked: set[str] | None = None) -> list[dict]:
    locked = locked or set()
    rows = []
    for p in PARAM_ORDER:
        lo, hi = DEFAULT_BOUNDS[p]
        v = params[p]
        frac = (np.log(v) - np.log(lo)) / (np.log(hi) - np.log(lo))
        rows.append({
            "param": p, "value": v, "lo": lo, "hi": hi,
            "log_frac": frac,
            "locked": p in locked,
            "at_bound": bool(p not in locked and (frac < 0.02 or frac > 0.98)),
        })
    return rows


def fmt_params(p: dict) -> str:
    return "  ".join(f"{k}={p[k]:.5g}" for k in PARAM_ORDER)


def run_country(cc: str, targets_path: str, outdir: str, args) -> dict:
    targets = json.load(open(targets_path))
    os.makedirs(outdir, exist_ok=True)
    log = []

    def say(s=""):
        print(s, flush=True)
        log.append(s)

    say("=" * 78)
    say(f"{cc}   targets: {targets_path}")
    say("=" * 78)
    say("targets: " + "  ".join(f"{k}={targets[k]:.4g}" for k in MOMENT_ORDER))

    fixed = COUNTRY_FIXED.get(cc, {})
    derive = COUNTRY_DERIVE.get(cc, {})
    locked = set(fixed) | set(derive)
    if locked:
        say(f"restricted: fixed={list(fixed)}  derived={list(derive)}  "
            f"-> {len(PARAM_ORDER) - len(locked)} free parameters")

    if args.quick:
        sim_kwargs = dict(n_workers=10_000, n_years_keep=36, steps_per_quarter=3)
        de_kw = dict(maxiter=8, popsize=4)
        tt_kw = dict(n_sobol=32, n_local=4)
        oos_seed = 555
    else:
        sim_kwargs = dict(n_workers=50_000, n_years_keep=36)
        de_kw = dict(maxiter=args.de_maxiter, popsize=args.de_popsize)
        tt_kw = dict(n_sobol=args.n_sobol, n_local=args.n_local)
        oos_seed = 987654

    # ---- differential evolution ----
    say("\n[DE] differential evolution ...")
    t0 = time.time()
    p_de, res_de = estimate(targets, sim_kwargs=sim_kwargs, workers=args.workers,
                            x0=KMV_TABLE3_PARAMS, disp=args.verbose,
                            fixed=fixed, derive=derive, **de_kw)
    say(f"[DE] done in {time.time()-t0:.0f}s   in-sample f={res_de.fun:.5f}")
    say("[DE] " + fmt_params(p_de))

    # ---- TikTak ----
    say("\n[TikTak] Sobol + sequential local ...")
    t0 = time.time()
    p_tt, info_tt = tiktak(targets, sim_kwargs=sim_kwargs, workers=args.workers,
                           disp=args.verbose, fixed=fixed, derive=derive,
                           x_start=KMV_TABLE3_PARAMS, **tt_kw)
    say(f"[TikTak] done in {time.time()-t0:.0f}s   in-sample f={info_tt['f_best']:.5f}")
    say("[TikTak] " + fmt_params(p_tt))

    # ---- out-of-sample scoring ----
    s_de, m_de = oos_score(p_de, targets, sim_kwargs, oos_seed)
    s_tt, m_tt = oos_score(p_tt, targets, sim_kwargs, oos_seed)
    say(f"\n[OOS] fresh-seed objective  DE={s_de:.5f}   TikTak={s_tt:.5f}")

    winner = "DE" if s_de <= s_tt else "TikTak"
    p_best = p_de if winner == "DE" else p_tt
    say(f"[OOS] canonical estimate: {winner}")

    # ---- agreement on the persistent component ----
    say("\n[persistent component] DE vs TikTak  (free params compared; "
        "locked params report the value implied by the two solutions)")
    pers = {}
    for k in PERSISTENT:
        a, b = p_de[k], p_tt[k]
        reldiff = abs(a - b) / max(abs(a), abs(b), 1e-12)
        tag = "  [locked]" if k in locked else ""
        pers[k] = dict(de=a, tiktak=b, reldiff=reldiff, locked=k in locked)
        say(f"  {k:8s} DE={a:.5g}   TikTak={b:.5g}   rel.diff={reldiff:.1%}{tag}")
    free_pers = [v["reldiff"] for k, v in pers.items() if k not in locked]
    max_pers_diff = max(free_pers) if free_pers else 0.0
    agree = max_pers_diff < 0.25
    say(f"  -> {'AGREE' if agree else 'DIVERGE'} on the free persistent-component "
        f"parameters (max rel.diff {max_pers_diff:.1%}; threshold 25%)")
    # ergodic sd of the persistent component under each solution
    for lab, pp in (("DE", p_de), ("TikTak", p_tt), ("canonical", p_best)):
        v2 = theoretical_stationary_var(pp["lambda2"], pp["beta2"], pp["sigma2"])
        say(f"  ergodic sd(z2) [{lab}] = {np.sqrt(v2):.3f}  (KMV: "
            f"{np.sqrt(KMV_ERGODIC_VAR2):.3f})")

    # ---- bound proximity ----
    say("\n[bounds] canonical estimate vs search box (log-scale position)")
    bt = bound_table(p_best, locked=locked)
    for r in bt:
        if r["locked"]:
            flag = "  [locked - not searched]"
        elif r["at_bound"]:
            flag = "  <-- AT BOUND"
        else:
            flag = ""
        say(f"  {r['param']:8s} {r['value']:10.5g}  in [{r['lo']:g}, {r['hi']:g}]"
            f"  pos={r['log_frac']:.2f}{flag}")
    any_bound = any(r["at_bound"] for r in bt)
    if any_bound:
        say("  -> WARNING: a free parameter sits at a bound; widen DEFAULT_BOUNDS and re-run.")
    else:
        say("  -> OK: no free parameter at a bound.")

    # ---- fit quality (canonical) ----
    say("\n[fit] 8-moment table (canonical estimate, OOS panel)")
    m_best = m_de if winner == "DE" else m_tt
    say(f"  {'moment':22s} {'data':>10s} {'model':>10s} {'rel.err':>9s}")
    for k in MOMENT_ORDER:
        re_ = (m_best[k] - targets[k]) / targets[k]
        say(f"  {MOMENT_LABELS[k]:22s} {targets[k]:10.4g} {m_best[k]:10.4g} {re_:>8.1%}")
    say(f"  sum of squared rel. errors (OOS): {rel_objective(m_best, targets):.5f}")

    # ---- regenerate table + chain from the canonical estimate ----
    say("\n[outputs] regenerating table2 + Markov chain (match_var_log=True) ...")
    tbl_sim = dict(n_workers=50_000, n_years_keep=36)
    disc_kwargs = dict(match_var_log=True, sim_kwargs=dict(n_workers=50_000))
    df, disc = build_table(targets, p_best, tbl_sim, disc_kwargs=disc_kwargs)
    df.to_csv(os.path.join(outdir, "table2.csv"), index=False)
    with open(os.path.join(outdir, "table2.tex"), "w") as f:
        f.write(df.to_latex(index=False, float_format="%.3f",
                            caption=f"Earnings Process Estimation Fit ({cc}, GRID 2001-2016)",
                            label=f"tab:earnings_fit_{cc.lower()}"))
    export_chain(disc, outdir)
    diag = diagnose(p_best, disc)
    say("  chain: " + json.dumps(diag["combined"]))
    for comp in ("1", "2"):
        d = diag[comp]
        say(f"  component {comp}: chain_var={d['chain_var']:.4g} "
            f"theory_var={d['theory_var']:.4g} n_states={d['n_states']}")
    say("  " + df.to_string(index=False).replace("\n", "\n  "))

    save_params(p_de, os.path.join(outdir, "params_de.json"))
    save_params(p_tt, os.path.join(outdir, "params_tiktak.json"))
    save_params(p_best, os.path.join(outdir, "params.json"))

    erg_sd2 = float(np.sqrt(theoretical_stationary_var(
        p_best["lambda2"], p_best["beta2"], p_best["sigma2"])))
    result = {
        "country": cc,
        "targets_path": targets_path,
        "targets": targets,
        "restricted": {"fixed": list(fixed), "derived": list(derive)},
        "params_de": p_de,
        "params_tiktak": p_tt,
        "params_canonical": p_best,
        "winner": winner,
        "in_sample_f": {"de": float(res_de.fun), "tiktak": float(info_tt["f_best"])},
        "oos_f": {"de": s_de, "tiktak": s_tt},
        "persistent_component": pers,
        "persistent_agree": bool(agree),
        "ergodic_sd_z2": {"canonical": erg_sd2, "kmv": float(np.sqrt(KMV_ERGODIC_VAR2))},
        "bounds": bt,
        "any_at_bound": bool(any_bound),
        "moments_canonical_oos": m_best,
        "chain_diagnostics": diag,
    }
    with open(os.path.join(outdir, "close_step1_report.json"), "w") as f:
        json.dump(result, f, indent=2, default=float)
    with open(os.path.join(outdir, "close_step1_log.txt"), "w") as f:
        f.write("\n".join(log))
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--countries", nargs="+", default=list(COUNTRIES),
                    choices=list(COUNTRIES))
    ap.add_argument("--outroot", default="output")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument("--de-maxiter", type=int, default=100)
    ap.add_argument("--de-popsize", type=int, default=12)
    ap.add_argument("--n-sobol", type=int, default=512)
    ap.add_argument("--n-local", type=int, default=10)
    args = ap.parse_args()

    t_all = time.time()
    summary = {}
    for cc in args.countries:
        outdir = os.path.join(args.outroot, cc.lower())
        summary[cc] = run_country(cc, COUNTRIES[cc], outdir, args)

    print("\n" + "#" * 78)
    print("SUMMARY")
    print("#" * 78)
    for cc, r in summary.items():
        rs = r["restricted"]
        print(f"\n{cc}: canonical={r['winner']}  "
              f"OOS f: DE={r['oos_f']['de']:.4f} TikTak={r['oos_f']['tiktak']:.4f}")
        if rs["fixed"] or rs["derived"]:
            print(f"   restricted: fixed={rs['fixed']} derived={rs['derived']}")
        print(f"   persistent component (free params): "
              f"{'AGREE' if r['persistent_agree'] else 'DIVERGE'}")
        print(f"   ergodic sd(z2): {r['ergodic_sd_z2']['canonical']:.2f} "
              f"(KMV {r['ergodic_sd_z2']['kmv']:.2f})")
        print(f"   free param at a bound: "
              f"{'YES -> widen bounds' if r['any_at_bound'] else 'no'}")
        print(f"   canonical params: {fmt_params(r['params_canonical'])}")
    with open(os.path.join(args.outroot, "close_step1_summary.json"), "w") as f:
        json.dump(summary, f, indent=2, default=float)
    print(f"\ntotal wall time: {time.time()-t_all:.0f}s")


if __name__ == "__main__":
    main()
