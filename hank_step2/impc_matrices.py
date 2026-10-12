"""
Task 5 (DECISIONS.md Sec 17): extract the intertemporal-MPC sequence M_{t,0}
(response of consumption at t to an unanticipated income shock at date 0 --
what Auclert-Rognlie-Straub's Figure 1 plots) and the capital-gains-response
sequence mcap_{t,0} = dC_t/dr_0, for every VALID chain, so these -- not the
direct/indirect decomposition alone -- are the primary cross-country/cross-
process comparison object per Sec 17's verdict.

"Valid chains" = the 4 canonical country estimates (usa, fra, ger, pol) plus
the two Poland sub-windows (pol_pre2009, pol_post2009). Explicitly EXCLUDES
every nested-specs variant chain (scripts/nested_specs.py), especially
beta2=0 (V3): a near-random-walk component's grid can blow up under
rescale_variance / match_var_log (discretize.py's new guard, added
alongside this script) and produce an artifactual ergodic distribution that
death.py would feed newborns from -- not something to run a HANK model on.

Reuses each chain's ALREADY-CALIBRATED (beta_hi, dbeta, omega) from
results/swap_death_<cc>.json rather than re-running the joint least_squares
calibration (expensive, already paid for) -- just rebuilds the steady state
once at the known optimum and takes Jacobians. Much cheaper than swap_death.py.

Usage: python impc_matrices.py
"""
from __future__ import annotations
import os, sys, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
AR = os.path.normpath(os.path.join(HERE, "..", "..", "annual-review"))
RESULTS = os.path.join(HERE, "results")
os.chdir(AR); sys.path.insert(0, AR); sys.path.insert(0, HERE)

import sequence_jacobian as sj
import household
from household import make_grids, income
from ge_blocks import (production, real_ST_bonds, fiscal, capitalization,
                       ex_post_r, nkpc, mkt_clearing)
from death import make_hh_death

ZETA = 1 / 180
T_IMPC = 40

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

VALID_CHAINS = ["usa", "fra", "ger", "pol", "pol_pre2009", "pol_post2009"]
# Poland ergodic-variance profile (Sec 23): computed separately via
# --extra, since two of these (0.3, 0.9103) do NOT converge in swap_death.py
# -- the Jacobian is still informative at the (non-equilibrium) point the
# optimiser landed at, but must be clearly flagged as such.
PROFILE_CHAINS = ["pol_var2_0.3", "pol_var2_0.4", "pol_var2_0.7", "pol_var2_0.9103"]


def load_chain(cc):
    from scipy import linalg
    d = os.path.normpath(os.path.join(HERE, "..", "output", cc))
    Q = np.loadtxt(f"{d}/income_process_Q.txt"); z = np.loadtxt(f"{d}/income_process_zgrid.txt")
    P = linalg.expm(Q); P /= P.sum(1)[:, None]
    assert np.allclose(P.sum(1), 1) and (P >= -1e-12).all()
    pi = sj.utilities.discretize.stationary(P); e = np.exp(z); e = e / (e @ pi)
    household.Pi_e, household.pi_e, household.e_grid_short, household.n_e = P, pi, e, len(pi)
    return float(pi @ (np.log(e) - pi @ np.log(e)) ** 2)


def compute_one(cc):
    cal_path = os.path.join(RESULTS, f"swap_death_{cc}.json")
    cal = json.load(open(cal_path))["calibration"]
    v = load_chain(cc)
    ha_params = hh_calib | dict(beta_hi=cal["beta_hi"], dbeta=cal["dbeta"],
                               omega=cal["omega"], beta_ave=cal["beta_ave"])
    ss = model_ha.steady_state({**ha_params, **common_params})
    J = hh.jacobian(ss, inputs=["wN_aftertax", "r"], outputs=["C"], T=T_IMPC)
    M_col0 = J["C", "wN_aftertax"][:, 0]
    mcap_col0 = J["C", "r"][:, 0]
    return {
        "country": cc, "ergodic_var_log_e": v,
        "beta_ave": cal["beta_ave"], "dbeta": cal["dbeta"],
        "converged": bool(cal["converged"]),
        "M_impc": [float(x) for x in M_col0],
        "mcap_impc": [float(x) for x in mcap_col0],
        "A": float(ss["A"]), "C": float(ss["C"]),
    }


def run_group(chains, json_name, png_name, title_suffix=""):
    out = {}
    for cc in chains:
        print(f"[impc] {cc} ...", flush=True)
        out[cc] = compute_one(cc)
        flag = "" if out[cc]["converged"] else "  [NOT CONVERGED -- off-equilibrium point]"
        print(f"  M_impc[0:4] = {[round(x,4) for x in out[cc]['M_impc'][:4]]}  "
              f"mcap_impc[0] = {out[cc]['mcap_impc'][0]:.4f}{flag}", flush=True)

    with open(os.path.join(RESULTS, json_name), "w") as f:
        json.dump(out, f, indent=2, default=float)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for cc, d in out.items():
        q = np.arange(len(d["M_impc"]))
        ls = ":" if not d["converged"] else "-"
        lbl = cc.upper() + ("" if d["converged"] else " (not converged)")
        axes[0].plot(q, d["M_impc"], label=lbl, lw=1.4, ms=3, marker="o", ls=ls)
        axes[1].plot(q, d["mcap_impc"], label=lbl, lw=1.4, ms=3, marker="o", ls=ls)
    axes[0].set_title(f"Intertemporal MPC: dC_t / dZ_0{title_suffix}")
    axes[1].set_title(f"Capital-gains response: dC_t / dr_0{title_suffix}")
    for ax in axes:
        ax.set_xlabel("quarters since shock"); ax.legend(fontsize=8); ax.axhline(0, color="grey", lw=0.6)
    axes[0].set_ylabel("consumption response")
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS, png_name), dpi=150)
    print(f"\nSaved: {RESULTS}/{json_name}, {png_name}")
    return out


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", action="store_true",
                    help="compute the Poland ergodic-variance profile chains instead of the 6 canonical ones")
    args = ap.parse_args()
    if args.profile:
        run_group(PROFILE_CHAINS, "impc_profile.json", "impc_profile.png",
                  " (Poland ergodic-variance profile)")
        return

    out = run_group(VALID_CHAINS, "impc_matrices.json", "impc_comparison.png")
    return out


if __name__ == "__main__":
    main()
