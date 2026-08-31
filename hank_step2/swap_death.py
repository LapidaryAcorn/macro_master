"""
Step 2 S5-redo: swap our chain into the one-asset model **with KMV stochastic
death added on top of ARS beta-heterogeneity** (which is what KMV actually have,
merged with what ARS actually have). Re-calibrate (beta_hi, dbeta, omega) jointly
to (A=20, impact labor MPC=0.20, SCF Lorenz). Report steady state, wealth / MPC
distribution, monetary + deficit IRFs, and the Fig-3b decomposition, all vs the
infinitely-lived baseline.

Usage:  python swap_death.py [USA|FRA|GER]   (default runs all three)

USA uses output/usa/ (the beta2-pinned chain, DECISIONS S5). TA/RA blocks do not
use the income process, so their IRFs are the baseline ones (not recomputed).
"""
from __future__ import annotations
import os, sys, json, time
import numpy as np
from scipy import optimize, linalg

HERE = os.path.dirname(os.path.abspath(__file__))
AR = os.path.normpath(os.path.join(HERE, "..", "..", "annual-review"))
os.chdir(AR); sys.path.insert(0, AR); sys.path.insert(0, HERE)

import sequence_jacobian as sj
import household
from household import make_grids, income
from ge_blocks import (production, real_ST_bonds, fiscal, capitalization,
                       ex_post_r, nkpc, mkt_clearing)
from death import make_hh_death

ZETA = 1 / 180
OUT = os.path.join(HERE, "results"); os.makedirs(OUT, exist_ok=True)
BASE = json.load(open(os.path.join(OUT, "baseline.json")))
lr = np.loadtxt("inputs/lorenz_nw_scf_2019.raw", delimiter=",")
pctl = np.arange(101) / 100
lorenz_scf = np.array([np.interp(p, lr[:, 0], lr[:, 1]) for p in pctl])

# GE calibration (chain-independent)
Y = 1.0; r = 0.005; A = 20.0; B = 4.0; G = 0.2; C = Y - G
tax = G + r * B; jf = A - B; div_post = r * jf
div = div_post / (1 - tax); w_post = (1 - tax) * (1 - div); mu = 1 / (1 - div)
eis = 1.0; frisch = 1.0
vscale = w_post / C ** (1 / eis) / Y ** (1 / frisch)
common_params = dict(Y=Y, r_ante=r, A=A, B=B, G=G, C=C, mu=mu, eis=eis, frisch=frisch,
                     vscale=vscale, pi=0, kappa=0.01, tax_rate_shock=0, T_shock=0,
                     zeta=0, T_rule_coeff=0)
hh_calib = dict(eis=1, min_a=0, max_a=4000, n_a=200, r=0.005, q=0.01, Tr_lumpsum=0,
                Tax_richest=0, zeta=0, N=1, wN_aftertax=0.7, zeta_d=ZETA, annuity=1.0)

hh = make_hh_death(ZETA).add_hetinputs([make_grids, income]); hh.name = "hh_ha"
common_blocks = [production, real_ST_bonds, fiscal, capitalization, ex_post_r, nkpc, mkt_clearing]
model_ha = sj.combine(common_blocks + [hh])


def load_chain(cc):
    d = os.path.normpath(os.path.join(HERE, "..", "output", cc))
    Q = np.loadtxt(f"{d}/income_process_Q.txt"); z = np.loadtxt(f"{d}/income_process_zgrid.txt")
    P = linalg.expm(Q); P /= P.sum(1)[:, None]
    assert np.allclose(P.sum(1), 1) and (P >= -1e-12).all()
    pi = sj.utilities.discretize.stationary(P); e = np.exp(z); e = e / (e @ pi)
    household.Pi_e, household.pi_e, household.e_grid_short, household.n_e = P, pi, e, len(pi)
    return float(pi @ (np.log(e) - pi @ np.log(e)) ** 2)


def lorenz(ss):
    D = ss.internals["hh_ha"]["D"].sum(0); ag = ss.internals["hh_ha"]["a_grid"]
    return np.array([np.interp(p, D.cumsum(), (ag * D).cumsum()) / ss["A"] for p in pctl]), D


def run(cc):
    t0 = time.time()
    v = load_chain(cc)
    R = {"country": cc.upper(), "with_death": True, "zeta": ZETA, "annuity": 1.0,
         "chain": dict(n_e=int(household.n_e), var_log_e_ergodic=v)}

    # ---- joint calibration (beta_hi, dbeta, omega) ----
    def resid(x):
        try:
            s = model_ha.steady_state({**(hh_calib | dict(beta_hi=x[0], dbeta=x[1], omega=x[2],
                                        beta_ave=x[0] - (1 - x[2]) * x[1])), **common_params})
        except ValueError:
            return np.array([1e2, 1e2, 1e2])
        m = hh.jacobian(s, inputs=["wN_aftertax"], outputs=["C"], T=2)["C", "wN_aftertax"][0, 0]
        lz, _ = lorenz(s)
        return np.array([s["A"] - 20, m - 0.20, (lz - lorenz_scf).sum()])

    sol = None
    for x0 in ([BASE["hh_params"]["beta_hi"], BASE["hh_params"]["dbeta"], BASE["hh_params"]["omega"]],
               [0.998, 0.05, 0.55], [0.9985, 0.12, 0.40]):
        s = optimize.least_squares(resid, x0, bounds=([0.90, 0.0, 0.02], [0.99999, 0.30, 0.98]),
                                   xtol=1e-11, ftol=1e-13)
        if sol is None or s.cost < sol.cost:
            sol = s
    bh, db, om = sol.x
    ba = bh - (1 - om) * db
    rr = resid(sol.x)
    b0 = BASE["hh_params"]
    R["calibration"] = dict(beta_hi=float(bh), beta_lo=float(bh - db), dbeta=float(db), omega=float(om),
                            beta_ave=float(ba), cost=float(sol.cost),
                            resid=dict(A=float(rr[0]), MPC=float(rr[1]), lorenz=float(rr[2])),
                            converged=bool(abs(rr[1]) < 0.01 and abs(rr[0]) < 0.3 and abs(rr[2]) < 0.3),
                            vs_baseline=dict(d_beta_ave=float(ba - b0["beta_ave"]),
                                             d_beta_lo=float((bh - db) - (b0["beta_hi"] - b0["dbeta"])),
                                             d_dbeta=float(db - b0["dbeta"])))

    # ---- steady state + wealth ----
    ha_params = hh_calib | dict(beta_hi=bh, dbeta=db, omega=om, beta_ave=ba)
    ss = model_ha.steady_state({**ha_params, **common_params})
    lz, D = lorenz(ss)
    R["steady_state"] = dict(A=float(ss["A"]), C=float(ss["C"]), asset_mkt=float(ss["asset_mkt"]),
                             goods_mkt=float(ss["goods_mkt"]))
    R["ha_wealth"] = dict(frac_at_borrowing_constraint=float(D[0]),
                          top10_share=float(1 - lz[90]), top1_share=float(1 - lz[99]),
                          lorenz_vs_scf_sum=float((lz - lorenz_scf).sum()),
                          lorenz=[float(x) for x in lz])

    # ---- MPCs ----
    T = 400
    Js = {"hh_ha": hh.jacobian(ss, inputs=["wN_aftertax", "N", "r"], outputs=["C", "A"], T=T)}
    M_lab = hh.jacobian(ss, inputs=["wN_aftertax"], outputs=["C"], T=26)["C", "wN_aftertax"]
    M_unw = hh.jacobian(ss, inputs=["Tr_lumpsum"], outputs=["C"], T=5)["C", "Tr_lumpsum"]
    R["mpc"] = dict(mpc_labor_impact=float(M_lab[0, 0]), mpc_unweighted_impact=float(M_unw[0, 0]),
                    share_spent_year1_unweighted=float((1 + r) ** (-np.arange(4)) @ M_unw[:4, 0]),
                    impc_labor_first8=[float(x) for x in M_lab[:8, 0]])

    # ---- monetary shock + decomposition ----
    dr = -0.25 * 0.9 ** np.arange(T)
    irf = model_ha.solve_impulse_linear(ss, unknowns=["Y", "B"],
            targets=["asset_mkt", "constant_owed_res"], inputs={"r_ante": dr},
            outputs=["Y", "r", "wN", "wN_aftertax"], Js=Js)
    R["monetary_irf_Y"] = [float(x) for x in irf["Y"][:41]]
    R["monetary_irf_Y_impact"] = float(irf["Y"][0])
    R["monetary_irf_Y_cum40"] = float(np.sum(irf["Y"][:40]))
    dC_cg = Js["hh_ha"]["C", "r"][:, 0] * irf["r"][0]
    dC_r = Js["hh_ha"]["C", "r"][:, 1:] @ irf["r"][1:]
    dC_l = Js["hh_ha"]["C", "wN_aftertax"] @ irf["wN"]
    dC_t = Js["hh_ha"]["C", "wN_aftertax"] @ (irf["wN_aftertax"] - irf["wN"])
    tot = irf["Y"]
    err = float(np.max(np.abs(dC_cg + dC_r + dC_l + dC_t - tot)))
    R["decomposition_check"] = dict(max_abs_error=err, rel=err / float(np.max(np.abs(tot))),
                                    ss_asset_mkt=float(ss["asset_mkt"]))

    def sh(n):
        T_ = float(np.sum(tot[:n]))
        return dict(direct_r=float(np.sum(dC_r[:n])) / T_,
                    indirect_labor=float(np.sum(dC_l[:n])) / T_,
                    indirect_tax=float(np.sum(dC_t[:n])) / T_,
                    indirect_cap_gains=float(np.sum(dC_cg[:n])) / T_,
                    indirect_total=float(np.sum(dC_l[:n] + dC_t[:n] + dC_cg[:n])) / T_)
    R["decomposition"] = dict(year1=sh(4), cum40=sh(40))

    # ---- deficit tax cut ----
    rho_B, rho = 0.975, 0.9
    dT = -rho ** np.arange(T); dBv = np.empty_like(dT); dBv[0] = -dT[0]
    for t in range(1, T):
        dBv[t] = rho_B * dBv[t - 1] - dT[t]
    irfB = model_ha.solve_impulse_linear(ss, unknowns=["Y"], targets=["asset_mkt"],
                                         inputs={"B": dBv}, outputs=["Y"], Js=Js)
    R["deficit_irf_Y_impact"] = float(irfB["Y"][0])
    R["deficit_irf_Y"] = [float(x) for x in irfB["Y"][:41]]

    # ---- vs baseline (infinitely-lived) ----
    def d(n, o): return dict(new=n, baseline=o, delta=n - o,
                             pct=(n - o) / o * 100 if o else float("nan"))
    R["vs_baseline"] = dict(
        beta_ave=d(ba, b0["beta_ave"]),
        frac_htm=d(R["ha_wealth"]["frac_at_borrowing_constraint"],
                   BASE["ha_wealth"]["frac_at_borrowing_constraint"]),
        top10=d(R["ha_wealth"]["top10_share"], BASE["ha_wealth"]["top10_share"]),
        mpc_labor_impact=d(R["mpc"]["mpc_labor_impact"], BASE["mpc"]["mpc_labor_impact"]),
        mpc_unweighted_impact=d(R["mpc"]["mpc_unweighted_impact"], BASE["mpc"]["mpc_unweighted_impact"]),
        monetary_Y_impact=d(R["monetary_irf_Y_impact"], BASE["monetary_irf_Y_impact"]["ha"]),
        monetary_Y_cum40=d(R["monetary_irf_Y_cum40"], BASE["monetary_irf_Y_cum40"]["ha"]),
        deficit_Y_impact=d(R["deficit_irf_Y_impact"], BASE["deficit_irf_Y_impact"]["ha"]),
        decomp_year1_direct=d(R["decomposition"]["year1"]["direct_r"],
                              BASE["decomposition"]["year1_shares"]["direct_r"]),
        decomp_year1_indirect=d(R["decomposition"]["year1"]["indirect_total"],
                                BASE["decomposition"]["year1_shares"]["indirect_total"]),
    )
    R["_secs"] = round(time.time() - t0)
    with open(os.path.join(OUT, f"swap_death_{cc}.json"), "w") as f:
        json.dump(R, f, indent=2)
    print(f"\n[{cc.upper()} + death] {R['_secs']}s  converged={R['calibration']['converged']}  "
          f"MPC={R['mpc']['mpc_labor_impact']:.3f}  frac_htm={R['ha_wealth']['frac_at_borrowing_constraint']:.3f}  "
          f"mon dY impact={R['monetary_irf_Y_impact']:.3f} (base {BASE['monetary_irf_Y_impact']['ha']:.3f})  "
          f"decomp yr1 {R['decomposition']['year1']['direct_r']*100:.0f}/{R['decomposition']['year1']['indirect_total']*100:.0f}",
          flush=True)
    return R


if __name__ == "__main__":
    ccs = [sys.argv[1].lower()] if len(sys.argv) > 1 else ["fra", "ger", "usa"]
    allR = {cc: run(cc) for cc in ccs}
    with open(os.path.join(OUT, "swap_death_summary.json"), "w") as f:
        json.dump(allR, f, indent=2)
