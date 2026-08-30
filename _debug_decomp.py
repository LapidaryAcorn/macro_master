"""Chase the ~10% decomposition residual with the death block.
Bisect: (1) goods market along the IRF path (dC == dY?), (2) Jacobian identity
(Js[C,r]@dr + Js[C,wNat]@dwNat == dC?), (3) shrink zeta_d toward 0 and watch the
residual -> if it stays ~10% the bug is structural in the overrides, not the
death magnitude."""
import os, sys, json, numpy as np
from scipy import optimize, linalg
AR = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "annual-review"))
os.chdir(AR); sys.path.insert(0, AR); sys.path.insert(0, os.path.dirname(__file__))
import sequence_jacobian as sj
import household
from household import make_grids, income, hh_ha
from ge_blocks import (production, real_ST_bonds, fiscal, capitalization,
                       ex_post_r, nkpc, mkt_clearing)
from death import make_hh_death

lr = np.loadtxt("inputs/lorenz_nw_scf_2019.raw", delimiter=","); pcts = np.arange(101)/100
lscf = np.array([np.interp(p, lr[:, 0], lr[:, 1]) for p in pcts])
Y = 1.0; r = 0.005; A = 20.0; B = 4.0; G = 0.2; C = Y - G
tax = G + r*B; jf = A - B; div_post = r*jf
div = div_post/(1-tax); w_post = (1-tax)*(1-div); mu = 1/(1-div)
vscale = w_post / C
common = dict(Y=Y, r_ante=r, A=A, B=B, G=G, C=C, mu=mu, eis=1.0, frisch=1.0, vscale=vscale,
              pi=0, kappa=0.01, tax_rate_shock=0, T_shock=0, zeta=0, T_rule_coeff=0)
blocks = [production, real_ST_bonds, fiscal, capitalization, ex_post_r, nkpc, mkt_clearing]


def load(cc):
    d = os.path.normpath(f"../kmv_grid_step1/output/{cc}")
    Q = np.loadtxt(f"{d}/income_process_Q.txt"); z = np.loadtxt(f"{d}/income_process_zgrid.txt")
    P = linalg.expm(Q); P /= P.sum(1)[:, None]
    pi = sj.utilities.discretize.stationary(P); e = np.exp(z); e = e/(e@pi)
    household.Pi_e, household.pi_e, household.e_grid_short, household.n_e = P, pi, e, len(pi)


def calib_and_decomp(hh, hc_extra, label, already_wired=False):
    hb = hh if already_wired else hh.add_hetinputs([make_grids, income])
    hb.name = "hh_ha"
    model = sj.combine(blocks + [hb])
    hc = dict(eis=1, min_a=0, max_a=4000, n_a=200, r=0.005, q=0.01, Tr_lumpsum=0,
              Tax_richest=0, zeta=0, N=1, wN_aftertax=0.7) | hc_extra

    def resid(x):
        try:
            s = model.steady_state({**(hc | dict(beta_hi=x[0], dbeta=x[1], omega=x[2],
                                     beta_ave=x[0]-(1-x[2])*x[1])), **common})
        except ValueError:
            return np.array([1e2, 1e2, 1e2])
        m = hb.jacobian(s, inputs=["wN_aftertax"], outputs=["C"], T=2)["C", "wN_aftertax"][0, 0]
        D = s.internals["hh_ha"]["D"].sum(0); ag = s.internals["hh_ha"]["a_grid"]
        lz = np.array([np.interp(p, D.cumsum(), (ag*D).cumsum())/s["A"] for p in pcts])
        return np.array([s["A"]-20, m-0.20, (lz-lscf).sum()])

    sol = optimize.least_squares(resid, [0.9985, 0.09, 0.55],
                                 bounds=([0.90, 0.0, 0.02], [0.99999, 0.30, 0.98]), xtol=1e-10)
    bh, db, om = sol.x
    ss = model.steady_state({**(hc | dict(beta_hi=bh, dbeta=db, omega=om,
                             beta_ave=bh-(1-om)*db)), **common})
    T = 300
    Js = {"hh_ha": hb.jacobian(ss, inputs=["wN_aftertax", "N", "r"], outputs=["C", "A"], T=T)}
    dr = -0.25 * 0.9 ** np.arange(T)
    irf = model.solve_impulse_linear(ss, unknowns=["Y", "B"],
            targets=["asset_mkt", "constant_owed_res"], inputs={"r_ante": dr},
            outputs=["Y", "C", "r", "wN", "wN_aftertax", "A"], Js=Js)
    dC_pred = Js["hh_ha"]["C", "r"] @ irf["r"] + Js["hh_ha"]["C", "wN_aftertax"] @ irf["wN_aftertax"]
    print(f"\n[{label}]  beta_lo={bh-db:.4f} dbeta={db:.4f}  resid={resid(sol.x)}")
    print(f"  goods mkt   max|dC - dY|          = {np.max(np.abs(irf['C']-irf['Y'])):.2e}   "
          f"(rel {np.max(np.abs(irf['C']-irf['Y']))/np.max(np.abs(irf['Y'])):.1%})")
    print(f"  Jac ident   max|Js@paths - dC|    = {np.max(np.abs(dC_pred-irf['C'])):.2e}   "
          f"(rel {np.max(np.abs(dC_pred-irf['C']))/np.max(np.abs(irf['C'])):.1%})")
    print(f"  asset resp  max|dA(hh) - dA(GE)|  = {np.max(np.abs(Js['hh_ha']['A','r']@irf['r'] + Js['hh_ha']['A','wN_aftertax']@irf['wN_aftertax'] - irf['A'])):.2e}")
    return irf, Js


if __name__ == "__main__":
    load("fra")
    calib_and_decomp(hh_ha, {}, "FRA no death", already_wired=True)
    for zd in [1/180, 1/1800, 1/18000]:
        calib_and_decomp(make_hh_death(zd), dict(zeta_d=zd, annuity=1.0), f"FRA death zeta_d=1/{round(1/zd)}")
    calib_and_decomp(make_hh_death(1/180), dict(zeta_d=1/180, annuity=0.0), "FRA death zeta_d=1/180 NO annuity")
