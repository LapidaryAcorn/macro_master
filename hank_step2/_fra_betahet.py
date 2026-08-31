"""Fast, focused: can a wider beta-heterogeneity get FRA (unrescaled chain) to
MPC 0.20 near A=20?  Point evaluations, no full calibration, no death."""
import os, sys, json, numpy as np
from scipy import linalg, optimize
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

def probe(cc, beta_hi, beta_lo, omega):
    ss = hh_ha.steady_state(hc | dict(beta_hi=beta_hi, dbeta=beta_hi - beta_lo, omega=omega))
    mpc = hh_ha.jacobian(ss, inputs=["wN_aftertax"], outputs=["C"], T=2)["C", "wN_aftertax"][0, 0]
    D = ss.internals["hh_ha"]["D"].sum(0); ag = ss.internals["hh_ha"]["a_grid"]
    lz = np.array([np.interp(p, D.cumsum(), (ag * D).cumsum()) / ss["A"] for p in pcts])
    return dict(A=float(ss["A"]), mpc=float(mpc), frac_a0=float(D[0]),
               lorenz_err=float((lz - lscf).sum()), top10=float(1 - lz[90]))

out = {}
for cc in ["fra", "ger"]:
    setup(cc)
    grid = []
    for beta_hi, beta_lo, omega in [
        (0.995, 0.90, 0.30), (0.995, 0.90, 0.50), (0.997, 0.85, 0.30),
        (0.997, 0.85, 0.50), (0.999, 0.80, 0.35), (0.999, 0.75, 0.40)]:
        try:
            r = probe(cc, beta_hi, beta_lo, omega)
        except Exception as ex:
            r = dict(error=str(type(ex).__name__))
        rec = dict(beta_hi=beta_hi, beta_lo=beta_lo, omega=omega, **r)
        grid.append(rec); print(cc, json.dumps(rec), flush=True)
    out[cc] = grid

json.dump(out, open(os.path.join(os.path.dirname(__file__), "results", "fra_betahet.json"), "w"), indent=2)
print("\nFor each country: does ANY (beta_lo, omega) reach MPC~0.20 with A near 20 and Lorenz near 0?")
for cc, g in out.items():
    good = [r for r in g if r.get("mpc", 0) > 0.15 and abs(r.get("A", 0) - 20) < 4]
    print(f"  {cc.upper()}: {'YES -> ' + json.dumps(good[0]) if good else 'no config in the probe reached MPC 0.15+ near A=20'}")
