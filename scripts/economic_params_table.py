"""
"Economically readable" version of the 6 raw quarterly-rate parameters, for
all 4 countries' canonical estimates: half-lives and mean inter-arrival
times in years (not quarterly rates), plus the ergodic sd each component
implies and the total. Saved to tables/10_economic_params.{csv,tex}.
"""
from __future__ import annotations
import json
import math
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# column label -> output dir. POL_V2 = Poland on the Sec 24 targets
# (var_log_earns from GRID std_loginc); POL = the v1 (std_perm) estimate.
COUNTRIES = {"USA": "usa", "FRA": "fra", "GER": "ger", "POL": "pol", "POL_V2": "pol_v2"}
ROWS = [
    ("Transitory half-life (years)", lambda p: math.log(2) / (4 * p["beta1"])),
    ("Persistent half-life (years)", lambda p: math.log(2) / (4 * p["beta2"])),
    ("Mean time between transitory shocks (years)", lambda p: 1 / (4 * p["lambda1"])),
    ("Mean time between persistent shocks (years)", lambda p: 1 / (4 * p["lambda2"])),
    ("sigma (transitory)", lambda p: p["sigma1"]),
    ("sigma (persistent)", lambda p: p["sigma2"]),
    ("Ergodic sd, transitory component",
     lambda p: math.sqrt(p["lambda1"] * p["sigma1"] ** 2 / (2 * p["beta1"]))),
    ("Ergodic sd, persistent component",
     lambda p: math.sqrt(p["lambda2"] * p["sigma2"] ** 2 / (2 * p["beta2"]))),
    ("Ergodic sd, combined (log earnings)",
     lambda p: math.sqrt(p["lambda1"] * p["sigma1"] ** 2 / (2 * p["beta1"])
                         + p["lambda2"] * p["sigma2"] ** 2 / (2 * p["beta2"]))),
]


def main():
    params = {lab: json.load(open(f"output/{d}/params.json")) for lab, d in COUNTRIES.items()}
    rows = []
    for label, fn in ROWS:
        rows.append({"quantity": label, **{cc: fn(p) for cc, p in params.items()}})
    df = pd.DataFrame(rows).set_index("quantity")

    os.makedirs("tables", exist_ok=True)
    df.round(4).to_csv("tables/10_economic_params.csv")
    with open("tables/10_economic_params.tex", "w") as f:
        f.write(df.to_latex(float_format="%.3f",
                            caption="Estimated earnings process: economically "
                                    "interpretable quantities (canonical estimates). "
                                    "USA and GER have beta2 pinned via ErgodicVarConstraint "
                                    "(Sec 8b); their persistent half-life/ergodic-sd rows "
                                    "reflect the calibration anchor, not a free estimate.",
                            label="tab:economic_params"))
    print(df.round(4).to_string())
    print("\nSaved: tables/10_economic_params.csv, tables/10_economic_params.tex")


if __name__ == "__main__":
    main()
