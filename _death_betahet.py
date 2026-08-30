"""Clean comparison: BOTH mechanisms active. Death (Blanchard-Yaari, annuities,
zeta=1/180) AND ARS beta-heterogeneity, with beta-het RE-CALIBRATED in the
presence of death. Full 3-target solve (A=20, impact labor MPC=0.20, SCF Lorenz),
multistart. Baseline chain + FRA + GER (USA is being re-estimated separately with
a pinned beta2, so its free-beta2 chain is moot)."""
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
hc = dict(eis=1, min_a=0, max_a=4000, n_a=200, r=0.005, q=0.01, Tr_lumpsum=0, Tax_richest=0,
          zeta=0, N=1, wN_aftertax=0.7, zeta_d=ZETA, annuity=1.0)
lr = np.loadtxt("inputs/lorenz_nw_scf_2019.raw", delimiter=","); pcts = np.arange(101)/100
lscf = np.array([np.interp(p, lr[:, 0], lr[:, 1]) for p in pcts])

def setup(cc):
    if cc == "baseline":
        P = linalg.expm(np.loadtxt("inputs/kmv_process/ymarkov_combined.txt")); P /= P.sum(1)[:, None]
        pi = sj.utilities.discretize.stationary(P)
        e = np.exp(np.loadtxt("inputs/kmv_process/ygrid_combined.txt")); e = e / (e @ pi)
    else:
        d = os.path.normpath(f"../kmv_grid_step1/output/{cc}")
        Q = np.loadtxt(f"{d}/income_process_Q.txt"); z = np.loadtxt(f"{d}/income_process_zgrid.txt")
        P = linalg.expm(Q); P /= P.sum(1)[:, None]
        pi = sj.utilities.discretize.stationary(P); e = np.exp(z); e = e / (e @ pi)
    household.Pi_e, household.pi_e, household.e_grid_short, household.n_e = P, pi, e, len(pi)

def resid(x):
    try:
        ss = hh.steady_state(hc | dict(beta_hi=x[0], dbeta=x[1], omega=x[2]))
    except ValueError:
        return np.array([1e2, 1e2, 1e2])
    mpc = hh.jacobian(ss, inputs=["wN_aftertax"], outputs=["C"], T=2)["C", "wN_aftertax"][0, 0]
    D = ss.internals["hh_ha"]["D"].sum(0); ag = ss.internals["hh_ha"]["a_grid"]
    lz = np.array([np.interp(p, D.cumsum(), (ag * D).cumsum()) / ss["A"] for p in pcts])
    return np.array([ss["A"] - 20, mpc - 0.20, (lz - lscf).sum()])

starts = [[0.9988, 0.088, 0.49], [0.999, 0.04, 0.55], [0.9995, 0.13, 0.35], [0.997, 0.06, 0.62]]
out = {}
for cc in ["baseline", "fra", "ger"]:
    setup(cc)
    best = None
    for x0 in starts:
        t = time.time()
        s = optimize.least_squares(resid, x0, bounds=([0.90, 0.0, 0.02], [0.99999, 0.40, 0.98]),
                                   xtol=1e-10, ftol=1e-12)
        r = resid(s.x)
        ss = hh.steady_state(hc | dict(beta_hi=s.x[0], dbeta=s.x[1], omega=s.x[2]))
        D = ss.internals["hh_ha"]["D"].sum(0)
        rec = dict(x0=x0, beta_hi=float(s.x[0]), beta_lo=float(s.x[0] - s.x[1]), dbeta=float(s.x[1]),
                   omega=float(s.x[2]), resid=[float(u) for u in r], cost=float(s.cost),
                   frac_a0=float(D[0]), A=float(ss["A"]), secs=round(time.time() - t))
        print(cc, json.dumps(rec), flush=True)
        if best is None or s.cost < best["cost"]:
            best = rec
    out[cc] = best

json.dump(out, open(os.path.join(os.path.dirname(__file__), "results", "death_betahet.json"), "w"), indent=2)
print("\n=== death + beta-heterogeneity, both re-calibrated (A=20, MPC=0.20, Lorenz) ===")
for cc, r in out.items():
    mpc = 0.20 + r["resid"][1]
    ok = abs(r["resid"][1]) < 0.02 and abs(r["resid"][0]) < 0.3 and abs(r["resid"][2]) < 0.3 and r["dbeta"] > 0.01
    print(f"{cc:9s}: {'CONVERGES' if ok else 'degenerate'}  MPC={mpc:.3f}  A={r['A']:.2f}  "
          f"dbeta={r['dbeta']:.4f} beta_lo={r['beta_lo']:.4f}  frac_a0={r['frac_a0']:.3f}  "
          f"lorenz_err={r['resid'][2]:+.2f}")
