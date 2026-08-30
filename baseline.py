"""
Step 2, task 1: reproduce the Auclert-Rognlie-Straub (2025, Annual Review)
one-asset HANK baseline UNCHANGED - their own KMV income process, their
calibration (hh_params.json) - and record the numbers everything else is judged
against.

Replicates the GE setup and the key results of `Annual Review main.ipynb`:
  - steady states of the HA / TA / RA models
  - household Jacobians
  - monetary policy shock IRFs (Fig 3a) for HA/TA/RA
  - the direct-vs-indirect decomposition (Fig 3b)
  - the deficit-financed tax cut IRFs (Fig 2a)

Run from anywhere; it chdirs into ../annual-review for the relative paths.
Writes results/baseline.json and results/baseline_*.png next to this file.
"""

from __future__ import annotations
import os, sys, json, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
AR = os.path.normpath(os.path.join(HERE, "..", "annual-review"))
os.chdir(AR)
sys.path.insert(0, AR)

import sequence_jacobian as sj
from household import hh_ha, hh_ta, hh_ra
import household

# ---------------------------------------------------------------------------
# GE blocks - copied verbatim from `Annual Review main.ipynb` section 1
# (they live in the notebook, not household.py).
# ---------------------------------------------------------------------------
@sj.simple
def production(Y, mu):
    N = Y
    wN = Y / mu
    div = Y - wN
    return N, wN, div

@sj.simple
def real_ST_bonds(r_ante):
    r_post_bonds = r_ante(-1)
    qbond = 1/(1 + r_ante)
    return r_post_bonds, qbond

@sj.simple
def nom_ST_bonds(r_ante, pi):
    pi_e = pi(1)
    r_post_bonds = (1 + r_ante(-1)) * (1 + pi_e(-1))/(1 + pi) - 1
    qbond = 1/(1 + r_ante)
    return r_post_bonds, qbond

@sj.simple
def fiscal(B, r_post_bonds, G, Y, wN, div, tax_rate_shock, T_shock, T_rule_coeff, qbond):
    T = (1 + r_post_bonds) * B(-1) + G - B
    tax_rate = T / Y
    wN_aftertax = (1 - tax_rate) * wN
    div_aftertax = (1 - tax_rate) * div
    tax_rate_res = tax_rate - tax_rate.ss - tax_rate_shock
    T_res = T - T.ss - T_rule_coeff * (B(-1) - B.ss) - T_shock
    constant_owed_res = B/qbond - B.ss/qbond.ss
    return T, tax_rate, wN_aftertax, div_aftertax, tax_rate_res, T_res, constant_owed_res

@sj.solved(unknowns={'p': (0.001,50)}, targets=['cap_cond'])
def capitalization(div_aftertax, r_ante, p):
    cap_cond = p - (div_aftertax(+1) + p(+1))/(1 + r_ante)
    r_post_equity = (div_aftertax + p)/p(-1) - 1
    return r_post_equity, cap_cond

@sj.simple
def ex_post_r(r_post_equity, r_post_bonds, B, p):
    r = (p(-1) * r_post_equity + B(-1) * r_post_bonds) / (p(-1) + B(-1))
    return r

@sj.simple
def nkpc(pi, kappa, N, Y, G, beta_ave, vscale, frisch, tax_rate, mu, eis):
    wedge = vscale*N**(1/frisch) - (1-tax_rate) / mu * (Y-G)**(-1/eis)
    nkpc_res = kappa * wedge + beta_ave * pi(+1) - pi
    return nkpc_res

@sj.simple
def taylor_rule(rstar, pi, phi):
    i = rstar + phi * pi
    r_ante = i - pi(+1)
    return r_ante

@sj.simple
def mkt_clearing(A, Y, C, p, B, G):
    asset_mkt = A - p - B
    goods_mkt = C + G - Y
    return asset_mkt, goods_mkt

OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)
R = {}   # everything we record

t0 = time.time()

# ---------------------------------------------------------------- income process
R["income_process"] = dict(
    source="their inputs/kmv_process/{ymarkov,ygrid}_combined.txt (KMV 2018)",
    n_e=int(household.n_e),
    Pi_e_rows_sum_to_1=bool(np.allclose(household.Pi_e.sum(1), 1)),
    mean_e=float(household.pi_e @ household.e_grid_short),
    var_log_e_ergodic=float(household.pi_e @ (np.log(household.e_grid_short)
                            - household.pi_e @ np.log(household.e_grid_short))**2),
)

# ---------------------------------------------------------------- GE calibration
# (verbatim from Annual Review main.ipynb, section 2)
Y = 1.0; r = 0.005; A = 20.0; B = 4.0; G = 0.2; C = Y - G
tax = G + r * B
j = A - B
div_post = r * j
div = div_post / (1 - tax)
w_post = (1 - tax) * (1 - div)
mu = 1 / (1 - div)
eis = 1.0; frisch = 1.0
vscale = w_post / C**(1/eis) / Y**(1/frisch)

common_params = dict(Y=Y, r_ante=r, A=A, B=B, G=G, C=C, mu=mu, eis=eis,
                     frisch=frisch, vscale=vscale, pi=0, kappa=0.01,
                     tax_rate_shock=0, T_shock=0, zeta=0, T_rule_coeff=0)
R["ge_calibration"] = dict(mu=float(mu), vscale=float(vscale), w_post=float(w_post),
                           div=float(div), tax=float(tax), r=r, A=A, B=B, G=G)

common_blocks = [production, real_ST_bonds, fiscal, capitalization, ex_post_r, nkpc, mkt_clearing]
models_hh = [hh_ha, hh_ta, hh_ra]
model_names = ["ha", "ta", "ra"]
models = {name: sj.combine(common_blocks + [hh]) for hh, name in zip(models_hh, model_names)}

# ---------------------------------------------------------------- steady states
ss = {}

ha_params = json.load(open("hh_params.json"))
R["hh_params"] = ha_params
ss["ha"] = models["ha"].steady_state({**ha_params, **common_params})

beta_ra = 1 / (1 + r)
ss["ra"] = models["ra"].steady_state({**dict(beta=beta_ra, beta_ave=beta_ra), **common_params},
                                     dissolve=["hh_ra"])
lam = 0.2 - r / (1 + r)
C_RA = (C - lam * w_post) / (1 - lam)
ss["ta"] = models["ta"].steady_state({**dict(beta=beta_ra, beta_ave=beta_ra, lam=lam, C_RA=C_RA),
                                      **common_params}, dissolve=["hh_ta"])

R["steady_state"] = {}
for k in model_names:
    R["steady_state"][k] = dict(
        A=float(ss[k]["A"]), C=float(ss[k]["C"]),
        asset_mkt=float(ss[k]["asset_mkt"]), goods_mkt=float(ss[k]["goods_mkt"]),
        nkpc_res=float(ss[k]["nkpc_res"]),
    )

# HA wealth distribution / Lorenz
D = ss["ha"].internals["hh_ha"]["D"].sum(axis=0)
a_grid = ss["ha"].internals["hh_ha"]["a_grid"]
pctl = np.arange(101) / 100
lorenz = np.array([np.interp(p, D.cumsum(), (a_grid * D).cumsum()) / ss["ha"]["A"] for p in pctl])
frac_hand_to_mouth = float(D[a_grid < 1e-9].sum()) if (a_grid < 1e-9).any() else float(D[0])
R["ha_wealth"] = dict(
    frac_at_borrowing_constraint=float(D[0]),
    lorenz_p50=float(lorenz[50]), lorenz_p90=float(lorenz[90]), lorenz_p99=float(lorenz[99]),
    top10_share=float(1 - lorenz[90]), top1_share=float(1 - lorenz[99]),
)

# ---------------------------------------------------------------- Jacobians
T = 400
Js = {"hh_ha": hh_ha.jacobian(ss["ha"], inputs=["wN_aftertax", "N", "r"],
                              outputs=["C", "A"], T=T)}

# intertemporal MPCs (Fig 1a) and aggregate MPC
M_labor = hh_ha.jacobian(ss["ha"], inputs=["wN_aftertax"], outputs=["C"], T=26)["C", "wN_aftertax"]
M_unwtd = hh_ha.jacobian(ss["ha"], inputs=["Tr_lumpsum"], outputs=["C"], T=5)["C", "Tr_lumpsum"]
R["mpc"] = dict(
    mpc_labor_impact=float(M_labor[0, 0]),
    mpc_unweighted_impact=float(M_unwtd[0, 0]),
    share_spent_year1_unweighted=float((1 + r) ** (-np.arange(4)) @ M_unwtd[:4, 0]),
    impc_labor_first8=[float(x) for x in M_labor[:8, 0]],
)

# ---------------------------------------------------------------- monetary shock
dr = -0.25 * 0.9 ** np.arange(T)
shock_r = {"r_ante": dr}
irfs_r = {k: models[k].solve_impulse_linear(
              ss[k], unknowns=["Y", "B"], targets=["asset_mkt", "constant_owed_res"],
              inputs=shock_r, outputs=["Y", "r", "wN", "wN_aftertax"], Js=Js)
          for k in model_names}
R["monetary_irf_Y"] = {k: [float(x) for x in irfs_r[k]["Y"][:41]] for k in model_names}
R["monetary_irf_Y_impact"] = {k: float(irfs_r[k]["Y"][0]) for k in model_names}
R["monetary_irf_Y_cum40"] = {k: float(np.sum(irfs_r[k]["Y"][:40])) for k in model_names}

# ---------------------------------------------------------------- decomposition (Fig 3b)
dC_cap_gains = Js["hh_ha"]["C", "r"][:, 0] * irfs_r["ha"]["r"][0]
dC_r = Js["hh_ha"]["C", "r"][:, 1:] @ irfs_r["ha"]["r"][1:]
dC_labor = Js["hh_ha"]["C", "wN_aftertax"] @ irfs_r["ha"]["wN"]
dC_tax = Js["hh_ha"]["C", "wN_aftertax"] @ (irfs_r["ha"]["wN_aftertax"] - irfs_r["ha"]["wN"])
total = irfs_r["ha"]["Y"]
assert np.allclose(dC_cap_gains + dC_r + dC_labor + dC_tax, total), "decomposition does not sum"

def _cum(x, n=40):
    return float(np.sum(x[:n]))

def _shares(n):
    tot = _cum(total, n)
    return dict(
        total=tot,
        direct_r=_cum(dC_r, n) / tot,
        indirect_labor=_cum(dC_labor, n) / tot,
        indirect_tax=_cum(dC_tax, n) / tot,
        indirect_cap_gains=_cum(dC_cap_gains, n) / tot,
        indirect_total=(_cum(dC_labor, n) + _cum(dC_tax, n) + _cum(dC_cap_gains, n)) / tot,
    )

R["decomposition"] = dict(
    note="monetary shock, HA model; Fig 3b of Annual Review. Shares comparable to "
         "KMV (2018) Table 7 col 1: direct r_b 19%, indirect w 51%, indirect T 32%, "
         "indirect r_a&q -2%  (=> 80% indirect / 20% direct, first-year average).",
    impact_levels=dict(total=float(total[0]), direct_r=float(dC_r[0]),
                       indirect_labor=float(dC_labor[0]), indirect_tax=float(dC_tax[0]),
                       indirect_cap_gains=float(dC_cap_gains[0])),
    impact_shares={k: (v / total[0] if k != "total" else float(total[0]))
                   for k, v in dict(total=total[0], direct_r=dC_r[0],
                                    indirect_labor=dC_labor[0], indirect_tax=dC_tax[0],
                                    indirect_cap_gains=dC_cap_gains[0]).items()},
    year1_shares=_shares(4),
    cum40_shares=_shares(40),
)

# ---------------------------------------------------------------- deficit tax cut (Fig 2a)
rho_B, rho = 0.975, 0.9
dT_shock = -rho ** np.arange(T)
dB = np.empty_like(dT_shock)
dB[0] = -dT_shock[0]
for t in range(1, T):
    dB[t] = rho_B * dB[t - 1] - dT_shock[t]
irfs_B = {k: models[k].solve_impulse_linear(ss[k], unknowns=["Y"], targets=["asset_mkt"],
                                            inputs={"B": dB}, outputs=["Y"], Js=Js)["Y"]
          for k in model_names}
R["deficit_irf_Y"] = {k: [float(x) for x in irfs_B[k][:41]] for k in model_names}
R["deficit_irf_Y_impact"] = {k: float(irfs_B[k][0]) for k in model_names}

R["_meta"] = dict(seconds=time.time() - t0, sj_from="github master",
                  python=sys.version.split()[0])

with open(os.path.join(OUT, "baseline.json"), "w") as f:
    json.dump(R, f, indent=2)

# ---------------------------------------------------------------- plots
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 3, figsize=(13, 3.4))
for k, c in zip(model_names, ["black", "orange", "red"]):
    ax[0].plot(irfs_r[k]["Y"][:41], label=k.upper(), color=c)
ax[0].set_title("Monetary shock: dY (Fig 3a)"); ax[0].legend(); ax[0].axhline(0, ls=":", c="gray")
for series, lab in [(total, "total"), (dC_r, "direct r"), (dC_labor, "indir. labor"),
                    (dC_tax, "indir. tax"), (dC_cap_gains, "indir. cap gains")]:
    ax[1].plot(series[:41], label=lab)
ax[1].set_title("Decomposition (Fig 3b)"); ax[1].legend(fontsize=7); ax[1].axhline(0, ls=":", c="gray")
for k, c in zip(model_names, ["black", "orange", "red"]):
    ax[2].plot(irfs_B[k][:41], label=k.upper(), color=c)
ax[2].set_title("Deficit tax cut: dY (Fig 2a)"); ax[2].legend(); ax[2].axhline(0, ls=":", c="gray")
plt.tight_layout()
plt.savefig(os.path.join(OUT, "baseline_irfs.png"), dpi=110)

print(json.dumps(R, indent=2))
print(f"\n[baseline] {R['_meta']['seconds']:.0f}s -> {OUT}/baseline.json")
