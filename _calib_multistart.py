"""Diagnostic: is the USA/FRA calibration failure a bad local min or structural?
Multi-start the (beta_hi, dbeta, omega) -> (A=20, MPC=0.20, SCF Lorenz) solve."""
import os, sys, json, numpy as np
from scipy import optimize, linalg
AR = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "annual-review"))
os.chdir(AR); sys.path.insert(0, AR)
import sequence_jacobian as sj
import household
from household import hh_ha

hc = dict(eis=1, min_a=0, max_a=4000, n_a=200, r=0.005, q=0.01, Tr_lumpsum=0,
          Tax_richest=0, zeta=0, N=1, wN_aftertax=0.7)
lr = np.loadtxt("inputs/lorenz_nw_scf_2019.raw", delimiter=","); pcts = np.arange(101) / 100
lscf = np.array([np.interp(p, lr[:, 0], lr[:, 1]) for p in pcts])

def setup(cc):
    d = os.path.normpath(f"../kmv_grid_step1/output/{cc}")
    Q = np.loadtxt(f"{d}/income_process_Q.txt"); z = np.loadtxt(f"{d}/income_process_zgrid.txt")
    P = linalg.expm(Q); P /= P.sum(1)[:, None]
    pi = sj.utilities.discretize.stationary(P)
    e = np.exp(z); e = e / (e @ pi)
    household.Pi_e, household.pi_e, household.e_grid_short, household.n_e = P, pi, e, len(pi)

def resid(x):
    try:
        ss = hh_ha.steady_state(hc | dict(beta_hi=x[0], dbeta=x[1], omega=x[2]))
    except ValueError:
        return np.array([1e3, 1e3, 1e3])
    mpc = hh_ha.jacobian(ss, inputs=['wN_aftertax'], outputs=['C'], T=2)['C', 'wN_aftertax'][0, 0]
    D = ss.internals['hh_ha']['D'].sum(0); ag = ss.internals['hh_ha']['a_grid']
    lz = np.array([np.interp(p, D.cumsum(), (ag * D).cumsum()) / ss['A'] for p in pcts])
    return np.array([ss['A'] - 20, mpc - 0.20, (lz - lscf).sum()])

starts = [[0.9988, 0.0879, 0.49], [0.995, 0.05, 0.5], [0.999, 0.12, 0.35],
          [0.9995, 0.20, 0.25], [0.99, 0.03, 0.65]]
out = {}
for cc in ["usa", "fra"]:
    setup(cc)
    best = None
    for x0 in starts:
        s = optimize.least_squares(resid, x0, bounds=([0.85, 0.0, 0.02], [0.9999, 0.35, 0.98]),
                                   xtol=1e-11, ftol=1e-13)
        r = resid(s.x)
        rec = dict(x0=x0, beta_hi=float(s.x[0]), dbeta=float(s.x[1]), omega=float(s.x[2]),
                   resid=[float(v) for v in r], cost=float(s.cost))
        print(cc, rec, flush=True)
        if best is None or s.cost < best["cost"]:
            best = rec
    out[cc] = best
json.dump(out, open(os.path.join(os.path.dirname(__file__), "results", "calib_multistart.json"), "w"), indent=2)
print("\nBEST:", json.dumps(out, indent=2))
