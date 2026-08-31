"""Is the USA/FRA S4 calibration failure structural, or an artifact of the
3-way least_squares collapsing dbeta? Fix dbeta on a grid; at each, calibrate
(beta_hi, omega) to (A=20, SCF Lorenz) and report the resulting MPC / frac-at-0.
No death. Unrescaled chains."""
import os, sys, json, numpy as np
from scipy import optimize, linalg
AR = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "annual-review"))
os.chdir(AR); sys.path.insert(0, AR); sys.path.insert(0, os.path.dirname(__file__))
import sequence_jacobian as sj
import household
from household import make_grids, income, hh_ha

lr = np.loadtxt("inputs/lorenz_nw_scf_2019.raw", delimiter=","); pcts = np.arange(101)/100
lscf = np.array([np.interp(p, lr[:, 0], lr[:, 1]) for p in pcts])
hc = dict(eis=1, min_a=0, max_a=4000, n_a=200, r=0.005, q=0.01, Tr_lumpsum=0, Tax_richest=0,
          zeta=0, N=1, wN_aftertax=0.7)

def setup(cc):
    d = os.path.normpath(f"../kmv_grid_step1/output/{cc}")
    Q = np.loadtxt(f"{d}/income_process_Q.txt"); z = np.loadtxt(f"{d}/income_process_zgrid.txt")
    P = linalg.expm(Q); P /= P.sum(1)[:, None]
    pi = sj.utilities.discretize.stationary(P); e = np.exp(z); e = e / (e @ pi)
    household.Pi_e, household.pi_e, household.e_grid_short, household.n_e = P, pi, e, len(pi)

out = {}
for cc in ["usa", "fra", "ger"]:
    setup(cc)
    rows = []
    for dbeta in [0.088, 0.13, 0.18, 0.25, 0.35]:
        def resid2(x):
            try:
                ss = hh_ha.steady_state(hc | dict(beta_hi=x[0], dbeta=dbeta, omega=x[1]))
            except ValueError:
                return np.array([1e3, 1e3])
            D = ss.internals["hh_ha"]["D"].sum(0); ag = ss.internals["hh_ha"]["a_grid"]
            lz = np.array([np.interp(p, D.cumsum(), (ag * D).cumsum()) / ss["A"] for p in pcts])
            return np.array([ss["A"] - 20, (lz - lscf).sum()])
        s = optimize.least_squares(resid2, [0.995, 0.5], bounds=([0.88, 0.02], [0.99999, 0.98]),
                                   xtol=1e-9, ftol=1e-11)
        ss = hh_ha.steady_state(hc | dict(beta_hi=s.x[0], dbeta=dbeta, omega=s.x[1]))
        mpc = hh_ha.jacobian(ss, inputs=["wN_aftertax"], outputs=["C"], T=2)["C", "wN_aftertax"][0, 0]
        D = ss.internals["hh_ha"]["D"].sum(0)
        r = resid2(s.x)
        rec = dict(dbeta=dbeta, beta_hi=float(s.x[0]), beta_lo=float(s.x[0] - dbeta),
                   omega=float(s.x[1]), A_err=float(r[0]), lorenz_err=float(r[1]),
                   mpc=float(mpc), frac_a0=float(D[0]))
        rows.append(rec)
        print(cc, json.dumps(rec), flush=True)
    out[cc] = rows

json.dump(out, open(os.path.join(os.path.dirname(__file__), "results", "betahet_sweep.json"), "w"), indent=2)
print("\n=== MPC achievable by widening beta-heterogeneity (A & Lorenz held) ===")
for cc, rows in out.items():
    best = max(rows, key=lambda r: r["mpc"] if abs(r["A_err"]) < 0.1 else -1)
    print(f"{cc.upper()}: max MPC = {best['mpc']:.3f} at dbeta={best['dbeta']} "
          f"(beta_lo={best['beta_lo']:.3f}), lorenz_err {best['lorenz_err']:+.2f}, frac_a0 {best['frac_a0']:.3f}")
