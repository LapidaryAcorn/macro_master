import os, sys, json, numpy as np
from scipy import optimize, linalg
AR = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "annual-review"))
os.chdir(AR); sys.path.insert(0, AR); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sequence_jacobian as sj
import household
from household import make_grids, income, hh_ha
lr = np.loadtxt("inputs/lorenz_nw_scf_2019.raw", delimiter=","); pcts = np.arange(101)/100
lscf = np.array([np.interp(p, lr[:, 0], lr[:, 1]) for p in pcts])
hc = dict(eis=1, min_a=0, max_a=4000, n_a=200, r=0.005, q=0.01, Tr_lumpsum=0, Tax_richest=0, zeta=0, N=1, wN_aftertax=0.7)
d = os.path.normpath("../kmv_grid_step1/output/fra")
Q = np.loadtxt(f"{d}/income_process_Q.txt"); z = np.loadtxt(f"{d}/income_process_zgrid.txt")
P = linalg.expm(Q); P /= P.sum(1)[:, None]; pi = sj.utilities.discretize.stationary(P); e = np.exp(z); e = e/(e@pi)
household.Pi_e, household.pi_e, household.e_grid_short, household.n_e = P, pi, e, len(pi)
def resid(x):
    try: ss = hh_ha.steady_state(hc | dict(beta_hi=x[0], dbeta=x[1], omega=x[2]))
    except ValueError: return np.array([1e2,1e2,1e2])
    mpc = hh_ha.jacobian(ss, inputs=['wN_aftertax'], outputs=['C'], T=2)['C','wN_aftertax'][0,0]
    D = ss.internals['hh_ha']['D'].sum(0); ag = ss.internals['hh_ha']['a_grid']
    lz = np.array([np.interp(p, D.cumsum(), (ag*D).cumsum())/ss['A'] for p in pcts])
    return np.array([ss['A']-20, mpc-0.20, (lz-lscf).sum()])
res = {}
for tag, x0 in [("from_GER_soln", [0.9943, 0.070, 0.50]), ("from_mid", [0.996, 0.05, 0.45]), ("from_wide", [0.998, 0.10, 0.40])]:
    s = optimize.least_squares(resid, x0, bounds=([0.90,0.005,0.05],[0.99998,0.30,0.95]), xtol=1e-10, ftol=1e-12)
    r = resid(s.x); ss = hh_ha.steady_state(hc | dict(beta_hi=s.x[0], dbeta=s.x[1], omega=s.x[2]))
    D = ss.internals['hh_ha']['D'].sum(0)
    rec = dict(x0=x0, beta_hi=float(s.x[0]), beta_lo=float(s.x[0]-s.x[1]), dbeta=float(s.x[1]), omega=float(s.x[2]),
               resid_A=float(r[0]), resid_MPC=float(r[1]), resid_lorenz=float(r[2]), cost=float(s.cost), frac_a0=float(D[0]))
    res[tag] = rec; print(json.dumps(rec), flush=True)
json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "fra_final.json"), "w"), indent=2)
ok = any(abs(r["resid_MPC"])<0.02 and abs(r["resid_A"])<0.3 and abs(r["resid_lorenz"])<0.3 and r["dbeta"]>0.01 for r in res.values())
print("\nFRA one-asset calibrates cleanly:", "YES" if ok else "NO - degenerate from every seed")
