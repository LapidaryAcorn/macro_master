"""
Extends the DECISIONS Sec 9 / output/tiktak_line_scan.json line scan PAST the
DE point, instead of only between TikTak and DE. Reuses the exact same
TikTak-recheck and DE points already on file (output/tiktak_line_scan.json),
the same log-parameter-space linear interpolation
(x(a) = exp((1-a) log(tiktak) + a log(de))), and the same objective (real
targets, default CRN seed, default sim_kwargs -- reproduced the original
a=0..1 points bit-for-bit as a check before extending).

a=0 is the TikTak point, a=1 is the DE (canonical) point, as in the original
file; this adds a=1.25, 1.5, 1.75, 2.0, i.e. continuing 25/50/75/100% of the
TikTak-to-DE distance *beyond* DE.

Question this answers: does the objective keep falling past the DE point in
the same direction, or does it turn around? If it keeps falling, the
canonical estimate is not even a local optimum along this direction and that
is a problem to know about now. If it turns around, the canonical point is
the valley floor along this direction -- the one-sentence result DECISIONS
Sec 9 wants for the write-up.

Usage: python scripts/extend_line_scan.py
"""
from __future__ import annotations

import json
import sys
import os

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from kmv_earnings.simulate import PARAM_ORDER
from kmv_earnings.estimate import make_objective
from kmv_earnings.grid_loader import load_targets

EXTEND_FRACTIONS = [1.25, 1.50, 1.75, 2.00]
TARGETS_PATH = {"FRA": "targets/fra_grid.json", "GER": "targets/ger_grid.json"}
ORIGINAL = json.load(open("output/tiktak_line_scan.json"))


def params_at(a: float, tiktak: dict, de: dict) -> dict:
    return {k: float(np.exp((1 - a) * np.log(tiktak[k]) + a * np.log(de[k])))
            for k in PARAM_ORDER}


def main():
    out = {}
    for cc, targets_path in TARGETS_PATH.items():
        targets = load_targets(targets_path)
        obj = make_objective(targets)
        tiktak = ORIGINAL[cc]["params_tiktak_recheck"]
        de = ORIGINAL[cc]["params_de"]

        # sanity check: reproduce a=0 and a=1 from the original file first
        for a_check, label in [(0.0, "tiktak"), (1.0, "de")]:
            p = params_at(a_check, tiktak, de)
            x = np.array([np.log(p[k]) for k in PARAM_ORDER])
            f = obj(x)
            f_orig = [pt["f"] for pt in ORIGINAL[cc]["line_scan"] if pt["a"] == a_check][0]
            print(f"[{cc}] reproduce a={a_check} ({label}): f={f:.6f}  "
                  f"(original {f_orig:.6f}, match={'YES' if abs(f-f_orig) < 1e-6 else 'NO -- investigate'})")

        extension = []
        for a in EXTEND_FRACTIONS:
            p = params_at(a, tiktak, de)
            x = np.array([np.log(p[k]) for k in PARAM_ORDER])
            f = obj(x)
            extension.append({"a": a, "beta1": p["beta1"], "params": p, "f": float(f)})
            print(f"[{cc}] a={a:.2f}  beta1={p['beta1']:.4f}  f={f:.6f}")

        f_de = [pt["f"] for pt in ORIGINAL[cc]["line_scan"] if pt["a"] == 1.0][0]
        f_vals = [pt["f"] for pt in extension]
        keeps_falling = all(f_vals[i] < f_vals[i - 1] for i in range(1, len(f_vals))) and f_vals[0] < f_de
        turns_around = f_vals[0] > f_de
        print(f"[{cc}] f at DE point (a=1.0) = {f_de:.6f}")
        print(f"[{cc}] verdict: {'KEEPS FALLING past DE -- canonical point is NOT a local optimum along this direction!' if keeps_falling else ('turns around / rises past DE -- canonical point is the valley floor along this direction' if turns_around else 'mixed / non-monotone past DE -- inspect')}")

        out[cc] = {
            "params_tiktak_recheck": tiktak, "params_de": de,
            "f_at_de": f_de, "extension": extension,
            "keeps_falling_past_de": bool(keeps_falling),
        }

    with open("output/tiktak_line_scan_extended.json", "w") as f:
        json.dump(out, f, indent=2, default=float)
    print("\nSaved: output/tiktak_line_scan_extended.json")


if __name__ == "__main__":
    main()
