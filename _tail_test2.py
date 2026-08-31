"""Test the USA tail hypothesis directly, cleanly: regenerate GER-shaped chains
from the KMV process family with the persistent-jump rate lambda2 scaled UP
(sigma2 scaled to hold the persistent component's ergodic variance) -> lower
1-year-change kurtosis, same variance. Then run the one-asset + death calibration.

  lambda2 up, sigma2 down (lam2*sig2^2 fixed): more frequent, smaller persistent
  jumps -> the annual change is more Gaussian -> lower kurtosis, unchanged var.
"""
import os, sys, json, numpy as np
from scipy import optimize, linalg
AR = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "annual-review"))
KMV = os.path.normpath(os.path.join(os.path.dirname(__file__), ".."))
os.chdir(AR); sys.path.insert(0, AR); sys.path.insert(0, os.path.dirname(__file__)); sys.path.insert(0, KMV)
import sequence_jacobian as sj
import household
from household import make_grids, income
from death import make_hh_death
from kmv_earnings.estimate import load_params
from kmv_earnings.discretize import discretize_process, simulate_discrete_panel
from kmv_earnings.simulate import moments_from_panel

ZETA = 1/180
lr = np.loadtxt("inputs/lorenz_nw_scf_2019.raw", delimiter=","); pcts = np.arange(101)/100
lscf = np.array([np.interp(p, lr[:, 0], lr[:, 1]) for p in pcts])
hh = make_hh_death(ZETA).add_hetinputs([make_grids, income]); hh.name = "hh_ha"
hc = dict(eis=1, min_a=0, max_a=4000, n_a=200, r=0.005, q=0.01, Tr_lumpsum=0, Tax_richest=0,
          zeta=0, N=1, wN_aftertax=0.7, zeta_d=ZETA, annuity=1.0)

p0 = load_params(os.path.join(KMV, "output/ger/params.json"))   # GER, calibrates cleanly

out = {}
for scale in [1.0, 2.0, 4.0, 8.0]:
    p = dict(p0)
    p["lambda2"] = p0["lambda2"] * scale
    p["sigma2"] = p0["sigma2"] / np.sqrt(scale)        # hold lambda2 * sigma2^2
    # keep beta2 so ergodic var of persistent comp lambda2 sigma2^2 / 2 beta2 is unchanged
    disc = discretize_process(p, match_var_log=True, sim_kwargs=dict(n_workers=40000))
    c = disc["combined"]
    P, pi, z = c["P"], c["pi"], c["zgrid"]
    m = moments_from_panel(simulate_discrete_panel(disc, n_workers=40000, n_years_keep=36))
    e = np.exp(z); e = e / (e @ pi)
    ergvar = float(pi @ (np.log(e) - pi @ np.log(e)) ** 2)
    household.Pi_e, household.pi_e, household.e_grid_short, household.n_e = P, pi, e, len(pi)

    def resid(x):
        try:
            s = hh.steady_state(hc | dict(beta_hi=x[0], dbeta=x[1], omega=x[2]))
        except ValueError:
            return np.array([1e2, 1e2, 1e2])
        mm = hh.jacobian(s, inputs=["wN_aftertax"], outputs=["C"], T=2)["C", "wN_aftertax"][0, 0]
        D = s.internals["hh_ha"]["D"].sum(0); ag = s.internals["hh_ha"]["a_grid"]
        lz = np.array([np.interp(pp, D.cumsum(), (ag*D).cumsum())/s["A"] for pp in pcts])
        return np.array([s["A"]-20, mm-0.20, (lz-lscf).sum()])

    b = optimize.least_squares(resid, [0.998, 0.067, 0.666],
                               bounds=([0.90, 0.0, 0.02], [0.99999, 0.30, 0.98]), xtol=1e-10)
    rr = resid(b.x)
    ss = hh.steady_state(hc | dict(beta_hi=b.x[0], dbeta=b.x[1], omega=b.x[2]))
    D = ss.internals["hh_ha"]["D"].sum(0)
    conv = abs(rr[1]) < 0.02 and abs(rr[0]) < 0.3 and abs(rr[2]) < 0.3 and b.x[1] > 0.01
    rec = dict(scale=scale, lambda2=p["lambda2"], sigma2=p["sigma2"],
               kurt_d1_disc=m["kurt_d1"], var_log_disc=m["var_log_earns"], ergodic_var_log_e=ergvar,
               dbeta=float(b.x[1]), beta_lo=float(b.x[0]-b.x[1]), mpc=0.20+float(rr[1]),
               frac_a0=float(D[0]), resid=[float(u) for u in rr], converged=bool(conv))
    out[scale] = rec
    print(json.dumps(rec), flush=True)

json.dump(out, open(os.path.join(os.path.dirname(__file__), "results", "tail_test2.json"), "w"), indent=2)
print("\n=== reduce kurtosis (hold variance) -> does the calibration break? ===")
for s, r in out.items():
    print(f"  lambda2 x{s}: kurt_d1(disc)={r['kurt_d1_disc']:5.2f}  ergvar={r['ergodic_var_log_e']:.2f}  "
          f"-> {'CONVERGES' if r['converged'] else 'COLLAPSES'}  (MPC {r['mpc']:.3f}, dbeta {r['dbeta']:.4f})")
