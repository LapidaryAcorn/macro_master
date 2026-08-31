"""Isolate the death-block Jacobian: does the LINEAR prediction Js @ shock_path
match the NONLINEAR impulse, for the household block alone (no GE)? Compare
no-death vs death vs death(zeta->0).  Fast: fixed calibration, one steady state
per case, T=150.
"""
import os, sys, json, numpy as np
from scipy import linalg
AR = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "annual-review"))
os.chdir(AR); sys.path.insert(0, AR); sys.path.insert(0, os.path.dirname(__file__))
import sequence_jacobian as sj
import household
from household import make_grids, income, hh_ha
from death import make_hh_death

# GER chain (calibrates cleanly)
d = os.path.normpath("../kmv_grid_step1/output/ger")
P = linalg.expm(np.loadtxt(f"{d}/income_process_Q.txt")); P /= P.sum(1)[:, None]
pi = sj.utilities.discretize.stationary(P)
e = np.exp(np.loadtxt(f"{d}/income_process_zgrid.txt")); e = e / (e @ pi)
household.Pi_e, household.pi_e, household.e_grid_short, household.n_e = P, pi, e, len(pi)

# GER + death calibration from _death_betahet.json
cal = dict(beta_hi=0.9980766933, dbeta=0.0665207932, omega=0.6660049212)
cal["beta_ave"] = cal["beta_hi"] - (1 - cal["omega"]) * cal["dbeta"]
base = dict(eis=1, min_a=0, max_a=4000, n_a=200, r=0.005, q=0.01, Tr_lumpsum=0, Tax_richest=0,
           zeta=0, N=1, wN_aftertax=0.7)
T = 150


def check(hb, extra, label, already=False):
    block = hb if already else hb.add_hetinputs([make_grids, income])
    block.name = "hh_ha"
    ss = block.steady_state({**base, **cal, **extra})
    J = block.jacobian(ss, inputs=["r", "wN_aftertax"], outputs=["C", "A"], T=T)
    out = {}
    for inp, amp in [("r", 1e-4), ("wN_aftertax", 1e-3)]:
        dpath = amp * 0.9 ** np.arange(T)
        lin = J["C", inp] @ dpath
        nl = block.impulse_nonlinear(ss, {inp: dpath}, outputs=["C"])["C"]
        rel = np.max(np.abs(lin - nl)) / np.max(np.abs(nl))
        out[inp] = dict(max_abs=float(np.max(np.abs(lin - nl))), rel=float(rel),
                        lin0=float(lin[0]), nl0=float(nl[0]))
    # also: A-response consistency (mass conservation under the shocked forward)
    dA_r = J["A", "r"] @ (1e-4 * 0.9 ** np.arange(T))
    nlA = block.impulse_nonlinear(ss, {"r": 1e-4 * 0.9 ** np.arange(T)}, outputs=["A"])["A"]
    out["A_r_rel"] = float(np.max(np.abs(dA_r - nlA)) / np.max(np.abs(nlA)))
    print(f"[{label}]  ss.A={ss['A']:.3f}  ss.C={ss['C']:.4f}")
    for k, v in out.items():
        if isinstance(v, dict):
            print(f"   Jac vs nonlinear, C wrt {k:12s}: rel {v['rel']:.2%}  (lin0 {v['lin0']:.3e} vs nl0 {v['nl0']:.3e})")
        else:
            print(f"   {k}: {v:.2%}")
    return out


R = {}
R["no_death"] = check(hh_ha, {}, "no death", already=True)
R["death_1_180"] = check(make_hh_death(1/180), dict(zeta_d=1/180, annuity=1.0), "death zeta=1/180")
R["death_tiny"] = check(make_hh_death(1e-6), dict(zeta_d=1e-6, annuity=1.0), "death zeta=1e-6 (should match no-death)")
R["death_noann"] = check(make_hh_death(1/180), dict(zeta_d=1/180, annuity=0.0), "death zeta=1/180 no annuity")
json.dump(R, open(os.path.join(os.path.dirname(__file__), "results", "decomp_jac_check.json"), "w"), indent=2)
