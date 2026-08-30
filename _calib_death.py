"""Does adding KMV stochastic death let the USA/FRA calibration converge on the
UNRESCALED chains? Re-run the (A=20, MPC=0.20, SCF Lorenz) solve with hh_ha_death."""
import os, sys, json, time, numpy as np
from scipy import optimize, linalg
AR = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "annual-review"))
os.chdir(AR); sys.path.insert(0, AR); sys.path.insert(0, os.path.dirname(__file__))
import sequence_jacobian as sj
import household
from household import make_grids, income
from death import make_hh_death

ZETA = 1/180
hh = make_hh_death(ZETA).add_hetinputs([make_grids, income]); hh.name = "hh_ha"
hc = dict(eis=1, min_a=0, max_a=4000, n_a=200, r=0.005, q=0.01, Tr_lumpsum=0,
          Tax_richest=0, zeta=0, N=1, wN_aftertax=0.7, zeta_d=ZETA)
lr = np.loadtxt("inputs/lorenz_nw_scf_2019.raw", delimiter=","); pcts = np.arange(101) / 100
lscf = np.array([np.interp(p, lr[:, 0], lr[:, 1]) for p in pcts])

def setup(cc):
    d = os.path.normpath(f"../kmv_grid_step1/output/{cc}")
    Q = np.loadtxt(f"{d}/income_process_Q.txt"); z = np.loadtxt(f"{d}/income_process_zgrid.txt")
    P = linalg.expm(Q); P /= P.sum(1)[:, None]
    pi = sj.utilities.discretize.stationary(P); e = np.exp(z); e = e / (e @ pi)
    household.Pi_e, household.pi_e, household.e_grid_short, household.n_e = P, pi, e, len(pi)
    return float(pi @ (np.log(e) - pi @ np.log(e)) ** 2)

def resid(x):
    try:
        ss = hh.steady_state(hc | dict(beta_hi=x[0], dbeta=x[1], omega=x[2]))
    except ValueError:
        return np.array([1e3, 1e3, 1e3])
    mpc = hh.jacobian(ss, inputs=['wN_aftertax'], outputs=['C'], T=2)['C', 'wN_aftertax'][0, 0]
    D = ss.internals['hh_ha']['D'].sum(0); ag = ss.internals['hh_ha']['a_grid']
    lz = np.array([np.interp(p, D.cumsum(), (ag * D).cumsum()) / ss['A'] for p in pcts])
    return np.array([ss['A'] - 20, mpc - 0.20, (lz - lscf).sum()])

starts = [[0.9988, 0.0879, 0.49], [0.9995, 0.05, 0.5], [0.9999, 0.15, 0.35], [0.999, 0.03, 0.65]]
out = {}
for cc in ["usa", "fra", "ger"]:
    v = setup(cc)
    best = None
    for x0 in starts:
        t = time.time()
        s = optimize.least_squares(resid, x0, bounds=([0.90, 0.0, 0.02], [0.99999, 0.4, 0.98]),
                                   xtol=1e-10, ftol=1e-12)
        r = resid(s.x)
        ss = hh.steady_state(hc | dict(beta_hi=s.x[0], dbeta=s.x[1], omega=s.x[2]))
        D = ss.internals['hh_ha']['D'].sum(0)
        rec = dict(x0=x0, beta_hi=float(s.x[0]), dbeta=float(s.x[1]), omega=float(s.x[2]),
                   beta_ave=float(s.x[0] - (1 - s.x[2]) * s.x[1]),
                   resid=[float(u) for u in r], cost=float(s.cost),
                   frac_a0=float(D[0]), A=float(ss['A']), secs=round(time.time() - t))
        print(cc, json.dumps(rec), flush=True)
        if best is None or s.cost < best["cost"]:
            best = rec
    best["chain_ergodic_var_log_e"] = v
    out[cc] = best

json.dump(out, open(os.path.join(os.path.dirname(__file__), "results", "calib_death.json"), "w"), indent=2)
print("\n=== BEST (with death, unrescaled chains) ===")
for cc, r in out.items():
    ok = abs(r["resid"][1]) < 0.01 and r["dbeta"] > 0.01
    print(f"{cc.upper()}: {'CONVERGES' if ok else 'still fails'}  "
          f"MPC={0.20 + r['resid'][1]:.3f}  A={r['A']:.2f}  dbeta={r['dbeta']:.4f}  "
          f"frac_a0={r['frac_a0']:.3f}  beta_ave={r['beta_ave']:.4f}  ergvar={r['chain_ergodic_var_log_e']:.2f}")
