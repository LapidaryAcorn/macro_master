"""Test the USA tail hypothesis DIRECTLY (not by elimination).

Take GER's chain (calibrates cleanly in one-asset + death). Compress its grid
tails with z -> sign(z)|z|^gamma, rescaled to hold var(log e) fixed, tuning
gamma so the 1-year-change kurtosis drops from GER's ~17.8 toward GRID-USA's
12.8. Then re-run the one-asset + death calibration.

  - collapses  -> the tail explanation for USA is confirmed; becomes a result.
  - still fine -> something else blocks USA; keep looking.
"""
import os, sys, json, numpy as np
from scipy import optimize, linalg
AR = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "annual-review"))
os.chdir(AR); sys.path.insert(0, AR); sys.path.insert(0, os.path.dirname(__file__))
import sequence_jacobian as sj
import household
from household import make_grids, income
from ge_blocks import (production, real_ST_bonds, fiscal, capitalization,
                       ex_post_r, nkpc, mkt_clearing)
from death import make_hh_death

ZETA = 1/180
lr = np.loadtxt("inputs/lorenz_nw_scf_2019.raw", delimiter=","); pcts = np.arange(101)/100
lscf = np.array([np.interp(p, lr[:, 0], lr[:, 1]) for p in pcts])
Y = 1.0; r = 0.005; A = 20.0; B = 4.0; G = 0.2; C = Y - G
tax = G + r*B; jf = A - B
div = (r*jf)/(1-tax); w_post = (1-tax)*(1-div); mu = 1/(1-div)
common = dict(Y=Y, r_ante=r, A=A, B=B, G=G, C=C, mu=mu, eis=1.0, frisch=1.0,
              vscale=w_post/C, pi=0, kappa=0.01, tax_rate_shock=0, T_shock=0, zeta=0, T_rule_coeff=0)
hh = make_hh_death(ZETA).add_hetinputs([make_grids, income]); hh.name = "hh_ha"
hc = dict(eis=1, min_a=0, max_a=4000, n_a=200, r=0.005, q=0.01, Tr_lumpsum=0, Tax_richest=0,
          zeta=0, N=1, wN_aftertax=0.7, zeta_d=ZETA, annuity=1.0)

d = os.path.normpath("../kmv_grid_step1/output/ger")
z0 = np.loadtxt(f"{d}/income_process_zgrid.txt")
P = linalg.expm(np.loadtxt(f"{d}/income_process_Q.txt")); P /= P.sum(1)[:, None]
pi = sj.utilities.discretize.stationary(P)
z0 = z0 - pi @ z0
var_target = float(pi @ z0**2)


def sim_kurt_d1(zg, P, pi, n=40000, years=36, seed=0):
    rng = np.random.default_rng(seed)
    cum = np.cumsum(P, axis=1)
    st = rng.choice(len(pi), size=n, p=pi)
    logy = np.empty((n, years)); af = np.zeros(n); yi = 0
    for q in range(4*years):
        u = rng.random(n); st = (u[:, None] > cum[st]).sum(1)
        af += np.exp(zg[st])
        if q % 4 == 3:
            logy[:, yi] = np.log(af); af[:] = 0; yi += 1
    d1 = (logy[:, 1:] - logy[:, :-1]).ravel(); d1 = d1 - d1.mean()
    return float(np.mean(d1**4) / np.mean(d1**2)**2), float(np.var(logy, axis=0).mean())


def compress(gamma):
    zc = np.sign(z0) * np.abs(z0) ** gamma
    zc *= np.sqrt(var_target / (pi @ (zc - pi @ zc) ** 2))
    return zc


base_k, base_v = sim_kurt_d1(z0, P, pi)
print(f"GER chain: kurt_d1(sim) = {base_k:.2f}   var_log = {base_v:.3f}   var_target(ergodic) = {var_target:.3f}", flush=True)

results = {}
for gamma in [1.0, 0.8, 0.65, 0.5, 0.42]:
    zc = compress(gamma)
    k, v = sim_kurt_d1(zc, P, pi)
    e = np.exp(zc); e = e / (e @ pi)
    household.Pi_e, household.pi_e, household.e_grid_short, household.n_e = P, pi, e, len(pi)

    def resid(x):
        try:
            s = hh.steady_state(hc | dict(beta_hi=x[0], dbeta=x[1], omega=x[2]))
        except ValueError:
            return np.array([1e2, 1e2, 1e2])
        m = hh.jacobian(s, inputs=["wN_aftertax"], outputs=["C"], T=2)["C", "wN_aftertax"][0, 0]
        D = s.internals["hh_ha"]["D"].sum(0); ag = s.internals["hh_ha"]["a_grid"]
        lz = np.array([np.interp(p, D.cumsum(), (ag*D).cumsum())/s["A"] for p in pcts])
        return np.array([s["A"]-20, m-0.20, (lz-lscf).sum()])

    best = optimize.least_squares(resid, [0.998, 0.067, 0.666],
                                  bounds=([0.90, 0.0, 0.02], [0.99999, 0.30, 0.98]), xtol=1e-10)
    rr = resid(best.x)
    ss = hh.steady_state(hc | dict(beta_hi=best.x[0], dbeta=best.x[1], omega=best.x[2]))
    D = ss.internals["hh_ha"]["D"].sum(0)
    conv = abs(rr[1]) < 0.02 and abs(rr[0]) < 0.3 and abs(rr[2]) < 0.3 and best.x[1] > 0.01
    rec = dict(gamma=gamma, kurt_d1=k, var_log_sim=v, dbeta=float(best.x[1]),
               beta_lo=float(best.x[0]-best.x[1]), mpc=0.20+float(rr[1]), frac_a0=float(D[0]),
               resid=[float(u) for u in rr], converged=bool(conv))
    results[gamma] = rec
    print(json.dumps(rec), flush=True)

json.dump(results, open(os.path.join(os.path.dirname(__file__), "results", "tail_test.json"), "w"), indent=2)
print("\n=== does reducing kurtosis (holding variance) break the calibration? ===")
for g, r_ in results.items():
    print(f"  gamma={g}: kurt_d1={r_['kurt_d1']:5.2f}  -> {'CONVERGES' if r_['converged'] else 'COLLAPSES'}  "
          f"(MPC {r_['mpc']:.3f}, dbeta {r_['dbeta']:.4f}, frac_a0 {r_['frac_a0']:.3f})")
