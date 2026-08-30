"""
Step 2, task S4/S3: swap our estimated Markov chain into the ARS one-asset HANK,
re-calibrate beta per chain, and compare to the baseline.

Convention (agreed): re-calibrate the discount-factor grid (beta_hi, dbeta,
omega) to the SAME targets ARS use - aggregate assets A=20, impact labor MPC
0.20, SCF Lorenz curve - so wealth matches by construction and the only thing
that can move steady-state MPCs / IRFs is the income process itself.

Usage:  python swap_chain.py [USA|FRA|GER]   (default USA)

Writes results/swap_<country>.json and results/swap_<country>.png; the JSON
carries a `vs_baseline` section with every delta.
"""

from __future__ import annotations
import os, sys, json, time
import numpy as np
from scipy import optimize, linalg

COUNTRY = (sys.argv[1] if len(sys.argv) > 1 else "USA").lower()

HERE = os.path.dirname(os.path.abspath(__file__))
AR = os.path.normpath(os.path.join(HERE, "..", "annual-review"))
CHAINDIR = os.path.normpath(os.path.join(HERE, "..", "kmv_grid_step1", "output", COUNTRY))
os.chdir(AR)
sys.path.insert(0, AR)

import sequence_jacobian as sj
import household
from household import hh_ha, hh_ta, hh_ra
from ge_blocks import (production, real_ST_bonds, fiscal, capitalization,
                       ex_post_r, nkpc, taylor_rule, mkt_clearing)

OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)
BASE = json.load(open(os.path.join(OUT, "baseline.json")))

# ---------------------------------------------------------------- swap the chain
Q = np.loadtxt(os.path.join(CHAINDIR, "income_process_Q.txt"))
zgrid = np.loadtxt(os.path.join(CHAINDIR, "income_process_zgrid.txt"))
Pi_new = linalg.expm(Q)
Pi_new /= Pi_new.sum(axis=1)[:, None]
assert np.allclose(Pi_new.sum(1), 1) and (Pi_new >= -1e-12).all(), "Pi not stochastic"
pi_new = sj.utilities.discretize.stationary(Pi_new)
e_new = np.exp(zgrid)
e_new = e_new / (e_new @ pi_new)

BASE_N_E = int(household.n_e)
BASE_VAR_LOG_E = float(household.pi_e @ (np.log(household.e_grid_short)
                       - household.pi_e @ np.log(household.e_grid_short)) ** 2)

# monkeypatch the module globals that make_grids / income read
household.Pi_e = Pi_new
household.pi_e = pi_new
household.e_grid_short = e_new
household.n_e = len(pi_new)

R = {"country": COUNTRY.upper(),
     "chain": dict(source=os.path.relpath(CHAINDIR, HERE), n_e=int(household.n_e),
                   var_log_e_ergodic=float(pi_new @ (np.log(e_new) - pi_new @ np.log(e_new)) ** 2),
                   mean_e=float(pi_new @ e_new))}

# ---------------------------------------------------------------- GE calibration (unchanged)
Y = 1.0; r = 0.005; A = 20.0; B = 4.0; G = 0.2; C = Y - G
tax = G + r * B; jf = A - B; div_post = r * jf
div = div_post / (1 - tax); w_post = (1 - tax) * (1 - div)
mu = 1 / (1 - div)
eis = 1.0; frisch = 1.0
vscale = w_post / C**(1/eis) / Y**(1/frisch)
common_params = dict(Y=Y, r_ante=r, A=A, B=B, G=G, C=C, mu=mu, eis=eis, frisch=frisch,
                     vscale=vscale, pi=0, kappa=0.01, tax_rate_shock=0, T_shock=0,
                     zeta=0, T_rule_coeff=0)
hh_calib = dict(eis=1, min_a=0, max_a=4000, n_a=200, r=0.005, q=0.01, Tr_lumpsum=0,
                Tax_richest=0, zeta=0, N=1, wN_aftertax=0.7)

# ---------------------------------------------------------------- re-calibrate beta
A_target, mpc_target = 20.0, 0.20
lorenz_raw = np.loadtxt("inputs/lorenz_nw_scf_2019.raw", delimiter=",")
percentiles = np.arange(101) / 100
lorenz_scf = np.array([np.interp(p, lorenz_raw[:, 0], lorenz_raw[:, 1]) for p in percentiles])

def get_lorenz(ss, pcts):
    D = ss.internals["hh_ha"]["D"].sum(axis=0)
    a_grid = ss.internals["hh_ha"]["a_grid"]
    return np.array([np.interp(p, D.cumsum(), (a_grid * D).cumsum()) / ss["A"] for p in pcts])

_neval = [0]
def resid(x):
    _neval[0] += 1
    try:
        ss = hh_ha.steady_state(hh_calib | dict(beta_hi=x[0], dbeta=x[1], omega=x[2]))
    except ValueError:                       # forward iteration didn't converge -> beta too high
        return np.array([1e3, 1e3, 1e3])
    mpc = hh_ha.jacobian(ss, inputs=["wN_aftertax"], outputs=["C"], T=2)["C", "wN_aftertax"][0, 0]
    lz = get_lorenz(ss, percentiles)
    return np.array([ss["A"] - A_target, mpc - mpc_target, (lz - lorenz_scf).sum()])

t0 = time.time()
b0_guess = [BASE["hh_params"]["beta_hi"], BASE["hh_params"]["dbeta"], BASE["hh_params"]["omega"]]
# bounded least-squares: beta_hi structurally < 1; our higher-variance chains need a lower beta
sol = optimize.least_squares(resid, b0_guess, bounds=([0.85, 0.0, 0.02], [0.9995, 0.35, 0.98]),
                             xtol=1e-10, ftol=1e-12, gtol=1e-12)
beta_hi, dbeta, omega = sol.x
beta_ave = beta_hi - (1 - omega) * dbeta
calib_secs = time.time() - t0
calib_resid = resid(sol.x)

b0 = BASE["hh_params"]
R["beta_recalibration"] = dict(
    n_evals=_neval[0], seconds=calib_secs, root_success=bool(sol.success),
    residuals_at_solution=dict(A_minus_target=float(calib_resid[0]),
                               mpc_minus_target=float(calib_resid[1]),
                               lorenz_sum_error=float(calib_resid[2])),
    new=dict(beta_hi=beta_hi, beta_lo=beta_hi - dbeta, dbeta=dbeta, omega=omega, beta_ave=beta_ave),
    baseline=dict(beta_hi=b0["beta_hi"], beta_lo=b0["beta_hi"] - b0["dbeta"],
                  dbeta=b0["dbeta"], omega=b0["omega"], beta_ave=b0["beta_ave"]),
    shift=dict(
        d_beta_hi=beta_hi - b0["beta_hi"],
        d_beta_lo=(beta_hi - dbeta) - (b0["beta_hi"] - b0["dbeta"]),
        d_beta_ave=beta_ave - b0["beta_ave"],
        d_omega=omega - b0["omega"],
        # annualised impatience of the average type, for intuition
        d_annual_discount_rate_ave=((1/beta_ave)**4 - 1) - ((1/b0["beta_ave"])**4 - 1),
    ),
)

# ---------------------------------------------------------------- steady states + GE
common_blocks = [production, real_ST_bonds, fiscal, capitalization, ex_post_r, nkpc, mkt_clearing]
model_ha = sj.combine(common_blocks + [hh_ha])
model_ta = sj.combine(common_blocks + [hh_ta])
model_ra = sj.combine(common_blocks + [hh_ra])

ha_params = dict(hh_calib) | dict(beta_hi=beta_hi, dbeta=dbeta, omega=omega, beta_ave=beta_ave)
ss_ha = model_ha.steady_state({**ha_params, **common_params})
beta_ra = 1 / (1 + r)
ss_ra = model_ra.steady_state({**dict(beta=beta_ra, beta_ave=beta_ra), **common_params}, dissolve=["hh_ra"])
lam = 0.2 - r / (1 + r); C_RA = (C - lam * w_post) / (1 - lam)
ss_ta = model_ta.steady_state({**dict(beta=beta_ra, beta_ave=beta_ra, lam=lam, C_RA=C_RA),
                               **common_params}, dissolve=["hh_ta"])
ss = dict(ha=ss_ha, ta=ss_ta, ra=ss_ra)

R["steady_state"] = {k: dict(A=float(ss[k]["A"]), C=float(ss[k]["C"]),
                             asset_mkt=float(ss[k]["asset_mkt"]),
                             goods_mkt=float(ss[k]["goods_mkt"]),
                             nkpc_res=float(ss[k]["nkpc_res"])) for k in ss}

D = ss_ha.internals["hh_ha"]["D"].sum(axis=0)
a_grid = ss_ha.internals["hh_ha"]["a_grid"]
pctl = np.arange(101) / 100
lz = np.array([np.interp(p, D.cumsum(), (a_grid * D).cumsum()) / ss_ha["A"] for p in pctl])
R["ha_wealth"] = dict(frac_at_borrowing_constraint=float(D[0]),
                      lorenz_p50=float(lz[50]), lorenz_p90=float(lz[90]), lorenz_p99=float(lz[99]),
                      top10_share=float(1 - lz[90]), top1_share=float(1 - lz[99]),
                      lorenz_vs_scf_sum=float((lz - lorenz_scf).sum()))

# ---------------------------------------------------------------- Jacobians / MPCs
T = 400
Js = {"hh_ha": hh_ha.jacobian(ss_ha, inputs=["wN_aftertax", "N", "r"], outputs=["C", "A"], T=T)}
M_labor = hh_ha.jacobian(ss_ha, inputs=["wN_aftertax"], outputs=["C"], T=26)["C", "wN_aftertax"]
M_unwtd = hh_ha.jacobian(ss_ha, inputs=["Tr_lumpsum"], outputs=["C"], T=5)["C", "Tr_lumpsum"]
R["mpc"] = dict(mpc_labor_impact=float(M_labor[0, 0]),
                mpc_unweighted_impact=float(M_unwtd[0, 0]),
                share_spent_year1_unweighted=float((1 + r) ** (-np.arange(4)) @ M_unwtd[:4, 0]),
                impc_labor_first8=[float(x) for x in M_labor[:8, 0]])

# ---------------------------------------------------------------- monetary shock
dr = -0.25 * 0.9 ** np.arange(T)
irfs_r = {k: m.solve_impulse_linear(ss[k], unknowns=["Y", "B"],
              targets=["asset_mkt", "constant_owed_res"], inputs={"r_ante": dr},
              outputs=["Y", "r", "wN", "wN_aftertax"], Js=Js)
          for k, m in dict(ha=model_ha, ta=model_ta, ra=model_ra).items()}
R["monetary_irf_Y"] = {k: [float(x) for x in irfs_r[k]["Y"][:41]] for k in irfs_r}
R["monetary_irf_Y_impact"] = {k: float(irfs_r[k]["Y"][0]) for k in irfs_r}
R["monetary_irf_Y_cum40"] = {k: float(np.sum(irfs_r[k]["Y"][:40])) for k in irfs_r}

# decomposition (Fig 3b)
dC_cap_gains = Js["hh_ha"]["C", "r"][:, 0] * irfs_r["ha"]["r"][0]
dC_r = Js["hh_ha"]["C", "r"][:, 1:] @ irfs_r["ha"]["r"][1:]
dC_labor = Js["hh_ha"]["C", "wN_aftertax"] @ irfs_r["ha"]["wN"]
dC_tax = Js["hh_ha"]["C", "wN_aftertax"] @ (irfs_r["ha"]["wN_aftertax"] - irfs_r["ha"]["wN"])
total = irfs_r["ha"]["Y"]
_decomp_sum = dC_cap_gains + dC_r + dC_labor + dC_tax
_decomp_err = float(np.max(np.abs(_decomp_sum - total)))
R["decomposition_check"] = dict(max_abs_error=_decomp_err,
                                rel_error=_decomp_err / float(np.max(np.abs(total))),
                                ss_asset_mkt=float(ss_ha["asset_mkt"]))
if _decomp_err > 1e-4 * np.max(np.abs(total)):
    print(f"WARNING: decomposition residual {_decomp_err:.2e} (rel {_decomp_err/np.max(np.abs(total)):.1e})")

def _cum(x, n): return float(np.sum(x[:n]))
def _shares(n):
    tot = _cum(total, n)
    return dict(total=tot, direct_r=_cum(dC_r, n)/tot, indirect_labor=_cum(dC_labor, n)/tot,
                indirect_tax=_cum(dC_tax, n)/tot, indirect_cap_gains=_cum(dC_cap_gains, n)/tot,
                indirect_total=(_cum(dC_labor, n)+_cum(dC_tax, n)+_cum(dC_cap_gains, n))/tot)
R["decomposition"] = dict(
    impact_shares={k: (v/total[0] if k != "total" else float(total[0])) for k, v in
                   dict(total=total[0], direct_r=dC_r[0], indirect_labor=dC_labor[0],
                        indirect_tax=dC_tax[0], indirect_cap_gains=dC_cap_gains[0]).items()},
    year1_shares=_shares(4), cum40_shares=_shares(40))

# deficit tax cut (Fig 2a)
rho_B, rho = 0.975, 0.9
dT = -rho ** np.arange(T)
dB = np.empty_like(dT); dB[0] = -dT[0]
for t in range(1, T):
    dB[t] = rho_B * dB[t-1] - dT[t]
irfs_B = {k: m.solve_impulse_linear(ss[k], unknowns=["Y"], targets=["asset_mkt"],
              inputs={"B": dB}, outputs=["Y"], Js=Js)["Y"]
          for k, m in dict(ha=model_ha, ta=model_ta, ra=model_ra).items()}
R["deficit_irf_Y_impact"] = {k: float(irfs_B[k][0]) for k in irfs_B}
R["deficit_irf_Y"] = {k: [float(x) for x in irfs_B[k][:41]] for k in irfs_B}

# ---------------------------------------------------------------- vs baseline
def _d(new, old): return dict(new=new, baseline=old, delta=new - old,
                              pct=(new - old) / old * 100 if old else float("nan"))
R["vs_baseline"] = dict(
    income_process=dict(n_e=_d(household.n_e, BASE_N_E),
                        var_log_e_ergodic=_d(R["chain"]["var_log_e_ergodic"], BASE_VAR_LOG_E)),
    beta_ave=_d(beta_ave, b0["beta_ave"]),
    beta_lo=_d(beta_hi - dbeta, b0["beta_hi"] - b0["dbeta"]),
    ha_frac_htm=_d(R["ha_wealth"]["frac_at_borrowing_constraint"],
                   BASE["ha_wealth"]["frac_at_borrowing_constraint"]),
    ha_top10=_d(R["ha_wealth"]["top10_share"], BASE["ha_wealth"]["top10_share"]),
    ha_top1=_d(R["ha_wealth"]["top1_share"], BASE["ha_wealth"]["top1_share"]),
    mpc_labor_impact=_d(R["mpc"]["mpc_labor_impact"], BASE["mpc"]["mpc_labor_impact"]),
    mpc_unweighted_impact=_d(R["mpc"]["mpc_unweighted_impact"], BASE["mpc"]["mpc_unweighted_impact"]),
    monetary_Y_impact_ha=_d(R["monetary_irf_Y_impact"]["ha"], BASE["monetary_irf_Y_impact"]["ha"]),
    monetary_Y_cum40_ha=_d(R["monetary_irf_Y_cum40"]["ha"], BASE["monetary_irf_Y_cum40"]["ha"]),
    deficit_Y_impact_ha=_d(R["deficit_irf_Y_impact"]["ha"], BASE["deficit_irf_Y_impact"]["ha"]),
    decomp_year1_direct=_d(R["decomposition"]["year1_shares"]["direct_r"],
                           BASE["decomposition"]["year1_shares"]["direct_r"]),
    decomp_year1_indirect=_d(R["decomposition"]["year1_shares"]["indirect_total"],
                             BASE["decomposition"]["year1_shares"]["indirect_total"]),
)

R["_meta"] = dict(seconds=time.time() - t0, python=sys.version.split()[0])
with open(os.path.join(OUT, f"swap_{COUNTRY}.json"), "w") as f:
    json.dump(R, f, indent=2)

# ---------------------------------------------------------------- plot vs baseline
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
ax[0].plot(pctl, lz, label=f"{COUNTRY.upper()} chain", color="black")
ax[0].plot(pctl, np.array(BASE["ha_wealth"].get("lorenz_full", [])) if "lorenz_full" in BASE["ha_wealth"] else lorenz_scf,
           label="SCF 2019", ls="--", color="green")
ax[0].plot([0, 1], [0, 1], ls=":", color="gray"); ax[0].set_title("HA wealth Lorenz"); ax[0].legend(fontsize=8)
ax[1].plot(R["monetary_irf_Y"]["ha"], label=f"{COUNTRY.upper()} chain", color="black")
ax[1].plot(BASE["monetary_irf_Y"]["ha"], label="baseline (their KMV)", ls="--", color="red")
ax[1].axhline(0, ls=":", c="gray"); ax[1].set_title("Monetary shock dY (HA)"); ax[1].legend(fontsize=8)
ax[2].plot(R["mpc"]["impc_labor_first8"], marker="o", label=f"{COUNTRY.upper()} chain", color="black")
ax[2].plot(BASE["mpc"]["impc_labor_first8"], marker="s", label="baseline", ls="--", color="red")
ax[2].set_title("iMPC out of labor income"); ax[2].legend(fontsize=8)
plt.tight_layout()
plt.savefig(os.path.join(OUT, f"swap_{COUNTRY}.png"), dpi=110)

print(json.dumps(R, indent=2))
print(f"\n[swap {COUNTRY.upper()}] {R['_meta']['seconds']:.0f}s -> results/swap_{COUNTRY}.json")
