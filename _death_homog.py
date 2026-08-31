"""KMV-style: homogeneous beta + stochastic death + annuities. Calibrate the
single beta to A=20, report MPC / frac-at-zero / Lorenz for the baseline chain
and our USA/FRA/GER chains (UNRESCALED)."""
import os, sys, json, numpy as np
from scipy import optimize, linalg
AR = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "annual-review"))
os.chdir(AR); sys.path.insert(0, AR); sys.path.insert(0, os.path.dirname(__file__))
import sequence_jacobian as sj
import household
from household import make_grids, income
from death import make_hh_death

ZETA = 1/180
lr = np.loadtxt("inputs/lorenz_nw_scf_2019.raw", delimiter=","); pcts = np.arange(101)/100
lscf = np.array([np.interp(p, lr[:, 0], lr[:, 1]) for p in pcts])
hh = make_hh_death(ZETA).add_hetinputs([make_grids, income]); hh.name = "hh_ha"
hc = dict(eis=1, min_a=0, max_a=4000, n_a=200, r=0.005, q=0.01, Tr_lumpsum=0, Tax_richest=0,
          zeta=0, N=1, wN_aftertax=0.7, zeta_d=ZETA, annuity=1.0)

def setup(cc):
    if cc == "baseline":
        P = linalg.expm(np.loadtxt("inputs/kmv_process/ymarkov_combined.txt"))
        P /= P.sum(1)[:, None]
        pi = sj.utilities.discretize.stationary(P)
        e = np.exp(np.loadtxt("inputs/kmv_process/ygrid_combined.txt")); e = e / (e @ pi)
    else:
        d = os.path.normpath(f"../kmv_grid_step1/output/{cc}")
        Q = np.loadtxt(f"{d}/income_process_Q.txt"); z = np.loadtxt(f"{d}/income_process_zgrid.txt")
        P = linalg.expm(Q); P /= P.sum(1)[:, None]
        pi = sj.utilities.discretize.stationary(P); e = np.exp(z); e = e / (e @ pi)
    household.Pi_e, household.pi_e, household.e_grid_short, household.n_e = P, pi, e, len(pi)
    return float(pi @ (np.log(e) - pi @ np.log(e)) ** 2)

out = {}
for cc in ["baseline", "usa", "fra", "ger"]:
    v = setup(cc)
    def Afn(b):
        return hh.steady_state(hc | dict(beta_hi=b, dbeta=0.0, omega=0.5))["A"] - 20
    try:
        b = optimize.brentq(Afn, 0.85, 0.99999, xtol=1e-10)
    except ValueError as ex:
        print(f"{cc}: brentq failed ({ex}) - trying wider");
        b = optimize.brentq(Afn, 0.5, 0.99999, xtol=1e-10)
    ss = hh.steady_state(hc | dict(beta_hi=b, dbeta=0.0, omega=0.5))
    mpc = hh.jacobian(ss, inputs=["wN_aftertax"], outputs=["C"], T=2)["C", "wN_aftertax"][0, 0]
    mpc_unw = hh.jacobian(ss, inputs=["Tr_lumpsum"], outputs=["C"], T=2)["C", "Tr_lumpsum"][0, 0]
    D = ss.internals["hh_ha"]["D"].sum(0); ag = ss.internals["hh_ha"]["a_grid"]
    lz = np.array([np.interp(p, D.cumsum(), (ag * D).cumsum()) / ss["A"] for p in pcts])
    rec = dict(chain_ergodic_var_log_e=v, beta=float(b), A=float(ss["A"]), mpc_labor=float(mpc),
               mpc_unweighted=float(mpc_unw), frac_a0=float(D[0]), lorenz_err=float((lz - lscf).sum()),
               top10=float(1 - lz[90]), top1=float(1 - lz[99]))
    out[cc] = rec
    print(cc, json.dumps(rec), flush=True)

json.dump(out, open(os.path.join(os.path.dirname(__file__), "results", "death_homog.json"), "w"), indent=2)
print("\n=== homogeneous-beta + death + annuities (calibrated to A=20) ===")
print(f"{'chain':9s} {'beta':>8s} {'MPC_lab':>8s} {'MPC_unw':>8s} {'frac_a0':>8s} {'lorenz_e':>9s} {'top10':>6s} {'ergvar':>7s}")
for cc, r in out.items():
    print(f"{cc:9s} {r['beta']:8.5f} {r['mpc_labor']:8.3f} {r['mpc_unweighted']:8.3f} "
          f"{r['frac_a0']:8.3f} {r['lorenz_err']:+9.2f} {r['top10']:6.2f} {r['chain_ergodic_var_log_e']:7.2f}")
