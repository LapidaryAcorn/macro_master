"""
Sec 24: robustness of a non-converged swap_death calibration. swap_death.py
tries 3 starting points for (beta_hi, dbeta, omega); when the best of them
stops at an INTERIOR point with a non-zero MPC residual, that may be a local
minimum of the least-squares problem rather than proof that the chain cannot
hit (A=20, MPC=0.20, SCF Lorenz). Re-run only the calibration step from
additional starts spanning the box (incl. the converged Polish v1 / var2=0.7
solutions and high-dbeta / low-omega corners) and report every endpoint.

Usage: python swap_death_multistart.py pol_v2
Writes results/swap_death_<cc>_multistart.json.
"""
from __future__ import annotations
import json
import os
import sys
import time

import numpy as np
from scipy import optimize

import swap_death as sd

STARTS = [
    [0.99792, 0.1252, 0.691],   # Polish v1 converged solution (Sec 20)
    [0.9975, 0.107, 0.607],     # Polish var2=0.7 converged solution (Sec 23)
    [0.995, 0.25, 0.50],
    [0.999, 0.20, 0.30],
    [0.996, 0.08, 0.80],
    [0.9995, 0.28, 0.70],
]


def main(cc):
    sd.load_chain(cc)

    def resid(x):
        try:
            s = sd.model_ha.steady_state({**(sd.hh_calib | dict(
                beta_hi=x[0], dbeta=x[1], omega=x[2], beta_ave=x[0] - (1 - x[2]) * x[1])),
                **sd.common_params})
        except ValueError:
            return np.array([1e2, 1e2, 1e2])
        m = sd.hh.jacobian(s, inputs=["wN_aftertax"], outputs=["C"], T=2)["C", "wN_aftertax"][0, 0]
        lz, _ = sd.lorenz(s)
        return np.array([s["A"] - 20, m - 0.20, (lz - sd.lorenz_scf).sum()])

    out = []
    for x0 in STARTS:
        t0 = time.time()
        s = optimize.least_squares(resid, x0, bounds=([0.90, 0.0, 0.02], [0.99999, 0.30, 0.98]),
                                   xtol=1e-11, ftol=1e-13)
        rr = resid(s.x)
        row = dict(x0=x0, x=[float(v) for v in s.x], cost=float(s.cost),
                   resid=dict(A=float(rr[0]), MPC=float(rr[1]), lorenz=float(rr[2])),
                   converged=bool(abs(rr[1]) < 0.01 and abs(rr[0]) < 0.3 and abs(rr[2]) < 0.3),
                   secs=time.time() - t0)
        out.append(row)
        print(json.dumps(row), flush=True)
        with open(os.path.join(sd.OUT, f"swap_death_{cc}_multistart.json"), "w") as f:
            json.dump(out, f, indent=2)


if __name__ == "__main__":
    main(sys.argv[1].lower())
