"""
Sec 31.3-31.5 tables from the finalized reports:
  * decomposition: v3 (small budget) and one KMV element at a time, plus the
    full KMV method and the full-budget v3 / v3_reset references;
  * Poland vs KMV, three versions (our method, KMV method, published chain).
Writes output/pol_v3_kmvmethod/tables.json and prints markdown.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from kmv_earnings.discretize import theoretical_stationary_var

ROOT = "output/pol_v3_kmvmethod"
LABEL = {
    "v3": "v3 (small budget)", "reset": "A1 reset jumps", "euler16": "A2–A4 16 steps/yr, Euler",
    "ergodic": "B1–B3, C1–C2 ergodic cross-section", "gap4": "C3 4-yr gap",
    "kurtw": "D2 kurtosis weight 0.5", "box": "E4 KMV box", "nsim5000": "B4 nsim 5000",
    "dfls": "E1–E3 DFLS from MakeGuess",
    "reset+ergodic": "A1 + B/C reset in the ergodic cross-section",
}


def load(path):
    return json.load(open(path)) if os.path.exists(path) else None


def best_by_start(rep):
    return {r["name"]: r for r in rep["runs"]}


def row(rep, reset=None):
    p, d = rep["params"], rep["derived"]
    ch = rep.get("chain_ours_5x15", {}).get("ergodic_var", {})
    pos = rep.get("bounds_log_position", {})
    at = [k for k, v in pos.items() if v < 0.005 or v > 0.995]
    return {"var2": d["var2"], "erg_sd_z2": d["erg_sd_z2"], "hl2": d["half_life_2_yr"],
            "chain_total": ch.get("total"), "chain_z2": ch.get("z2"),
            "f_own": rep["oos_f"], "f_eq": rep.get("f_eq"), "at_bound": at,
            "canonical": rep.get("canonical"), "params": p}


def main():
    out = {"decomp": {}, "refs": {}}
    for v in LABEL:
        rep = load(os.path.join(ROOT, "decomp", v, "report.json"))
        if rep:
            out["decomp"][v] = row(rep)
            out["decomp"][v]["starts"] = {n: {"oos_f": r["oos_f"], "var2": r["derived"]["var2"]}
                                          for n, r in best_by_start(rep).items()}
    k = load(os.path.join(ROOT, "kmv", "report.json"))
    if k:
        out["kmv"] = row(k)
        out["kmv"]["chain_kmv_3x11"] = k.get("chain_kmv_3x11", {}).get("ergodic_var")
        out["kmv"]["starts"] = {n: {"oos_f": r["oos_f"], "var2": r["derived"]["var2"],
                                    "params": r["params"]} for n, r in best_by_start(k).items()}
    us = load(os.path.join(ROOT, "us_v3", "report.json"))
    if us:
        out["us_v3"] = row(us)
    for name, path, reset in (("v3 full", "output/pol_v3/report.json", False),
                              ("v3_reset full", "output/pol_v3_reset/report.json", True)):
        rep = load(path)
        if rep:
            p = rep["params"]
            out["refs"][name] = {"var2": theoretical_stationary_var(p["lambda2"], p["beta2"], p["sigma2"], reset),
                                 "oos_f": rep["oos_f"], "params": p}
    json.dump(out, open(os.path.join(ROOT, "tables.json"), "w"), indent=2, default=float)

    f = lambda x, n=3: "–" if x is None else f"{x:.{n}f}"
    print("| element | var₂ | ergodic sd z₂ | half-life z₂ (yr) | chain var(log e), 5×15 | f (own weights, OOS) | f (equal weights) | at bound |")
    print("|---|---|---|---|---|---|---|---|")
    for v, r in out["decomp"].items():
        print(f"| {LABEL[v]} | {f(r['var2'])} | {f(r['erg_sd_z2'])} | {f(r['hl2'],1)} | {f(r['chain_total'])} | "
              f"{f(r['f_own'],4)} | {f(r['f_eq'],4)} | {', '.join(r['at_bound']) or '–'} |")
    if k:
        r = out["kmv"]
        print(f"| **all KMV elements (§31.3)** | {f(r['var2'])} | {f(r['erg_sd_z2'])} | {f(r['hl2'],1)} | "
              f"{f(r['chain_total'])} | {f(r['f_own'],4)} | {f(r['f_eq'],4)} | {', '.join(r['at_bound']) or '–'} |")
    print()
    for n, r in out.get("kmv", {}).get("starts", {}).items():
        print(n, r)
    print(json.dumps(out.get("us_v3"), default=float, indent=1))


if __name__ == "__main__":
    main()
