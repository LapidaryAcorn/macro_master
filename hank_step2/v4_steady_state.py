"""
DECISIONS Sec 32 (item 8): one-asset ARS HANK with KMV stochastic death
(swap_death.py: ARS beta-heterogeneity + death zeta = 1/180, annuity), steady
state only, on the Poland v4 chains:
  output/pol_v4/estimate/chain_ours_5x15   (our 5 x 15 scheme)
  output/pol_v4/estimate/chain_kmv_3x11    (KMV's 3 x 11 procedure)
Joint calibration of (beta_hi, dbeta, omega) to A = 20, impact labour MPC 0.20
and the SCF Lorenz curve, exactly as swap_death.run (same starts, same bounds).
No Jacobians beyond the MPC, no IRFs.

Usage: python hank_step2/v4_steady_state.py
Output: hank_step2/results/v4_steady_state.json
"""
from __future__ import annotations
import json, os, sys, time
import numpy as np
from scipy import optimize

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import swap_death as sd          # chdir to annual-review + model set-up

CHAINS = {"v4_5x15": "pol_v4/estimate/chain_ours_5x15",
          "v4_3x11": "pol_v4/estimate/chain_kmv_3x11"}


def mpc(s):
    return sd.hh.jacobian(s, inputs=["wN_aftertax"], outputs=["C"], T=2)["C", "wN_aftertax"][0, 0]


def run(name, cc):
    t0 = time.time()
    v = sd.load_chain(cc)
    R = {"chain": cc, "n_e": int(sd.household.n_e), "var_log_e_ergodic": v,
         "zeta": sd.ZETA, "model": "one-asset ARS + KMV death (swap_death.py)"}
    n_ss = [0]

    def ss_of(x):
        n_ss[0] += 1
        return sd.model_ha.steady_state({**(sd.hh_calib | dict(beta_hi=x[0], dbeta=x[1], omega=x[2],
                                         beta_ave=x[0] - (1 - x[2]) * x[1])), **sd.common_params})

    def resid(x):
        try:
            s = ss_of(x)
        except ValueError:
            return np.array([1e2, 1e2, 1e2])
        lz, _ = sd.lorenz(s)
        return np.array([s["A"] - 20, mpc(s) - 0.20, (lz - sd.lorenz_scf).sum()])

    b = sd.BASE["hh_params"]
    sol, starts = None, []
    for x0 in ([b["beta_hi"], b["dbeta"], b["omega"]], [0.998, 0.05, 0.55], [0.9985, 0.12, 0.40]):
        s = optimize.least_squares(resid, x0, bounds=([0.90, 0.0, 0.02], [0.99999, 0.30, 0.98]),
                                   xtol=1e-11, ftol=1e-13)
        starts.append({"x0": list(x0), "x": s.x.tolist(), "cost": float(s.cost)})
        if sol is None or s.cost < sol.cost:
            sol = s
    t_cal = time.time() - t0
    bh, db, om = sol.x
    rr = resid(sol.x)
    t1 = time.time()
    ss = ss_of(sol.x)
    t_ss = time.time() - t1
    lz, D = sd.lorenz(ss)
    ag = ss.internals["hh_ha"]["a_grid"]
    cdf = np.cumsum(D)
    q = {f"p{p}": float(np.interp(p / 100, cdf, ag)) for p in (10, 25, 50, 75, 90, 99)}
    mean_a = float(ag @ D)
    R.update(
        calibration=dict(beta_hi=float(bh), beta_lo=float(bh - db), dbeta=float(db), omega=float(om),
                         beta_ave=float(bh - (1 - om) * db),
                         beta_hi_annual=float(bh ** 4), beta_lo_annual=float((bh - db) ** 4),
                         cost=float(sol.cost), resid=dict(A=float(rr[0]), MPC=float(rr[1]), lorenz=float(rr[2])),
                         converged=bool(abs(rr[1]) < 0.01 and abs(rr[0]) < 0.3 and abs(rr[2]) < 0.3),
                         starts=starts),
        mpc_labor_impact=float(mpc(ss)),
        wealth=dict(A=float(ss["A"]), mean=mean_a, frac_at_borrowing_constraint=float(D[0]),
                    frac_below_1q_income=float(np.interp(sd.hh_calib["wN_aftertax"], ag, cdf)),
                    top10_share=float(1 - lz[90]), top1_share=float(1 - lz[99]),
                    bottom50_share=float(lz[50]), quantiles=q,
                    gini=float(1 - 2 * np.trapezoid(lz, sd.pctl)),
                    lorenz_vs_scf_sum=float((lz - sd.lorenz_scf).sum()), lorenz=[float(x) for x in lz]),
        market_clearing=dict(asset_mkt=float(ss["asset_mkt"]), goods_mkt=float(ss["goods_mkt"])),
        timing=dict(calibration_s=round(t_cal, 1), n_steady_states=n_ss[0],
                    one_steady_state_s=round(t_ss, 2), total_s=round(time.time() - t0, 1)))
    print(f"[{name}] var(log e)={v:.3f} beta_ave={R['calibration']['beta_ave']:.5f} "
          f"conv={R['calibration']['converged']} MPC={R['mpc_labor_impact']:.3f} "
          f"HtM={D[0]:.3f} top10={R['wealth']['top10_share']:.3f} t={R['timing']['total_s']}s", flush=True)
    return R


if __name__ == "__main__":
    out = {"baseline_infinitely_lived": {k: sd.BASE[k] for k in ("hh_params", "ha_wealth") if k in sd.BASE}}
    for name, cc in CHAINS.items():
        out[name] = run(name, cc)
    path = os.path.join(sd.OUT, "v4_steady_state.json")
    json.dump(out, open(path, "w"), indent=2)
    print("wrote", path)
