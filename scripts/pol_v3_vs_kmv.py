"""
Sec 28 (E, F): Poland v3 vs the KMV process and chain; chain grid-size test.

E. Table: parameters; economic quantities (half-lives, years between shocks,
   ergodic sds under the additive formula lambda sigma^2 / (2 beta) used in this
   repo AND under KMV's published reset reading lambda sigma^2 / (lambda + 2 beta),
   Sec 29.5); moments of the continuous process (ages 25-55, and 25-60 for
   KMV's own design); moments of the chains; transition-matrix properties of
   the annual chain P^4 (probability of staying in the same state, expected
   stay in the top / bottom decile of e, annual autocorrelation of log e);
   ergodic var(log e).
   Chains: Poland v3 (output/pol_v3, 5x15), KMV Table 3 re-discretized by us
   (5x15), and KMV's published chain (annual-review/inputs/kmv_process, 3x11).

F. Poland v3 chains at 3x11, 3x15, 5x11, 5x15: chain vs continuous moments
   (ages 25-55), and the time of one steady-state solve of the one-asset HA
   household block of hank_step2 (ARS hh_ha at the baseline beta parameters,
   no recalibration; median of 5 solves after one warm-up).

Usage: python scripts/pol_v3_vs_kmv.py {table|grid} [--params ...] [--chain ...]
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.linalg import expm

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from kmv_earnings.simulate import MOMENT_ORDER, KMV_TABLE3_PARAMS, model_moments
from kmv_earnings.discretize import (discretize_process, simulate_discrete_panel,
                                     theoretical_stationary_var)
from kmv_earnings.simulate import moments_from_panel

N_AGE = 31
AR = os.path.normpath(os.path.join(ROOT, "..", "annual-review"))
KMV_CHAIN = os.path.join(AR, "inputs", "kmv_process")
OUT = "output/pol_v3/vs_kmv"


def econ(p):
    r = {}
    for j in ("1", "2"):
        lam, b, s = p[f"lambda{j}"], p[f"beta{j}"], p[f"sigma{j}"]
        r[f"half_life_{j}_yr"] = np.log(2) / b / 4
        r[f"years_between_shocks_{j}"] = 1 / lam / 4
        r[f"erg_sd_{j}_additive"] = np.sqrt(theoretical_stationary_var(lam, b, s))
        r[f"erg_sd_{j}_reset"] = np.sqrt(lam * s ** 2 / (lam + 2 * b))
    r["erg_var_total_additive"] = r["erg_sd_1_additive"] ** 2 + r["erg_sd_2_additive"] ** 2
    r["erg_var_total_reset"] = r["erg_sd_1_reset"] ** 2 + r["erg_sd_2_reset"] ** 2
    return r


def chain_from_files(zgrid_path, q_path):
    zg = np.loadtxt(zgrid_path); Q = np.loadtxt(q_path)
    return zg, Q


def chain_props(zg, Q):
    P = expm(Q); P /= P.sum(1, keepdims=True)
    w, v = np.linalg.eig(P.T)
    pi = np.abs(np.real(v[:, np.argmin(np.abs(w - 1))])); pi /= pi.sum()
    le = np.log(np.exp(zg) / (pi @ np.exp(zg)))
    P4 = np.linalg.matrix_power(P, 4)
    m = pi @ le; var = pi @ (le - m) ** 2
    ac = (pi * (le - m)) @ P4 @ (le - m) / var
    order = np.argsort(le)

    def decile_stay(idx_sorted, top):
        # states making up (about) the top/bottom 10% of ergodic mass; a state
        # straddling the cut-off enters with its partial mass
        idx = idx_sorted[::-1] if top else idx_sorted
        cum = np.cumsum(pi[idx]); wgt = np.clip(0.10 - (cum - pi[idx]), 0, pi[idx])
        wgt = wgt / wgt.sum()
        stay_q = wgt @ (P[np.ix_(idx, idx)] @ (wgt > 0).astype(float))  # quarterly prob. to remain
        return float(1 / (1 - stay_q) / 4)                                # expected stay, years
    return {"n_states": len(zg), "ergodic_var_log_e": float(var),
            "prob_stay_same_state_annual": float(pi @ np.diag(P4)),
            "prob_stay_same_state_quarterly": float(pi @ np.diag(P)),
            "exp_years_top_decile": decile_stay(order, True),
            "exp_years_bottom_decile": decile_stay(order, False),
            "autocorr_log_e_annual": float(ac),
            "mass_on_extreme_states": float(pi[order[0]] + pi[order[-1]])}


def chain_moments(disc_or_none, zg=None, Q=None, n_workers=100_000):
    if disc_or_none is None:
        # wrap a combined chain as a single 'component' for the panel simulator
        P = expm(Q); P /= P.sum(1, keepdims=True)
        w, v = np.linalg.eig(P.T)
        pi = np.abs(np.real(v[:, np.argmin(np.abs(w - 1))])); pi /= pi.sum()
        z0 = int(np.argmin(np.abs(zg)))
        # z=0 start: the published chain has no exact zero state; use the state
        # nearest to the ergodic mean of log e instead
        z0 = int(np.argmin(np.abs(zg - pi @ zg)))
        disc = {"1": {"grid": zg, "P": P, "pi": pi},
                "2": {"grid": np.zeros(1), "P": np.ones((1, 1)), "pi": np.ones(1)}}
        logy = _sim_from_state(disc, z0, n_workers)
    else:
        logy = simulate_discrete_panel(disc_or_none, n_workers=n_workers, n_years_keep=36)
    return moments_from_panel(logy, n_age_cols=N_AGE)


def _sim_from_state(disc, z0, n_workers, seed=4321):
    rng = np.random.default_rng(seed)
    cum = np.cumsum(disc["1"]["P"], axis=1)
    s = np.full(n_workers, z0)
    logy = np.empty((n_workers, 36)); flow = np.zeros(n_workers)
    for q in range(4 * 36):
        u = rng.random(n_workers)
        s = (u[:, None] > cum[s]).sum(axis=1)
        flow += np.exp(disc["1"]["grid"][s])
        if q % 4 == 3:
            logy[:, q // 4] = np.log(flow); flow[:] = 0
    return logy


def cmd_table(a):
    os.makedirs(a.out, exist_ok=True)
    pl = json.load(open(a.params))
    kmv = KMV_TABLE3_PARAMS
    rows = {}
    for lab, p in (("Poland v3", pl), ("KMV Table 3", kmv)):
        rows[lab] = {"params": p, "econ": econ(p),
                     "cont_moments_25_55": model_moments(p, n_workers=100_000, n_years_keep=36,
                                                         n_age_cols=N_AGE, seed=987654),
                     "cont_moments_25_60": model_moments(p, n_workers=100_000, n_years_keep=36,
                                                         seed=987654)}
    zg, Q = chain_from_files(os.path.join(a.chain, "income_process_zgrid.txt"),
                             os.path.join(a.chain, "income_process_Q.txt"))
    d_pl = discretize_process(pl, match_var_log=True, sim_kwargs=dict(n_workers=50_000, n_age_cols=N_AGE))
    if not np.allclose(d_pl["combined"]["zgrid"], zg, atol=1e-6):
        print("WARNING: exported chain differs from the rebuilt one; using the exported grid")
    d_kmv = discretize_process(kmv, match_var_log=True, sim_kwargs=dict(n_workers=50_000, n_age_cols=N_AGE))
    zk, Qk = chain_from_files(os.path.join(KMV_CHAIN, "ygrid_combined.txt"),
                              os.path.join(KMV_CHAIN, "ymarkov_combined.txt"))
    chains = {"Poland v3 chain (5x15)": (d_pl, zg, Q),
              "KMV Table 3, our chain (5x15)": (d_kmv, d_kmv["combined"]["zgrid"], d_kmv["combined"]["Q"]),
              "KMV published chain (3x11)": (None, zk, Qk)}
    for lab, (d, z_, Q_) in chains.items():
        rows[lab] = {"chain_props": chain_props(z_, Q_),
                     "chain_moments_25_55": chain_moments(d, z_, Q_)}
    json.dump(rows, open(os.path.join(a.out, "vs_kmv.json"), "w"), indent=2, default=float)

    def line(name, vals, fmt="{:.4g}"):
        return f"| {name} | " + " | ".join(fmt.format(v) if isinstance(v, (int, float, np.floating)) else str(v) for v in vals) + " |"
    L = ["| | Poland v3 | KMV Table 3 |", "|---|---|---|"]
    for k in ("lambda1", "beta1", "sigma1", "lambda2", "beta2", "sigma2"):
        L.append(line(k, [rows["Poland v3"]["params"][k], kmv[k]]))
    for k in rows["Poland v3"]["econ"]:
        L.append(line(k, [rows["Poland v3"]["econ"][k], rows["KMV Table 3"]["econ"][k]]))
    for k in MOMENT_ORDER:
        L.append(line(f"continuous {k} (25-55)", [rows["Poland v3"]["cont_moments_25_55"][k],
                                                 rows["KMV Table 3"]["cont_moments_25_55"][k]]))
    for k in MOMENT_ORDER:
        L.append(line(f"continuous {k} (25-60)", [rows["Poland v3"]["cont_moments_25_60"][k],
                                                 rows["KMV Table 3"]["cont_moments_25_60"][k]]))
    L += ["", "| chain | " + " | ".join(chains) + " |", "|---|" + "---|" * len(chains)]
    for k in rows["Poland v3 chain (5x15)"]["chain_props"]:
        L.append(line(k, [rows[c]["chain_props"][k] for c in chains]))
    for k in MOMENT_ORDER:
        L.append(line(f"chain {k} (25-55)", [rows[c]["chain_moments_25_55"][k] for c in chains]))
    open(os.path.join(a.out, "vs_kmv.md"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("\n".join(L))


def ss_timer():
    """Return a function chain(zgrid, Q) -> seconds for one hh_ha steady state."""
    cwd = os.getcwd()
    os.chdir(AR); sys.path.insert(0, AR)
    import household
    from household import hh_ha
    base = json.load(open(os.path.join(ROOT, "hank_step2", "results", "baseline.json")))
    os.chdir(cwd)
    hh_calib = dict(eis=1, min_a=0, max_a=4000, n_a=200, r=0.005, q=0.01, Tr_lumpsum=0,
                    Tax_richest=0, zeta=0, N=1, wN_aftertax=0.7)
    hp = base["hh_params"]
    calib = hh_calib | dict(beta_hi=hp["beta_hi"], dbeta=hp["dbeta"], omega=hp["omega"])

    def run(zg, Q, reps=5):
        P = expm(Q); P /= P.sum(1, keepdims=True)
        w, v = np.linalg.eig(P.T)
        pi = np.abs(np.real(v[:, np.argmin(np.abs(w - 1))])); pi /= pi.sum()
        e = np.exp(zg); e /= e @ pi
        household.Pi_e, household.pi_e, household.e_grid_short, household.n_e = P, pi, e, len(pi)
        try:
            hh_ha.steady_state(calib)                       # warm-up (numba compile etc.)
            ts = []
            for _ in range(reps):
                t0 = time.perf_counter(); ss = hh_ha.steady_state(calib); ts.append(time.perf_counter() - t0)
            return float(np.median(ts)), float(ss["A"]), None
        except Exception as ex:                             # e.g. no convergence at baseline betas
            return None, None, repr(ex)[:200]
    return run


def cmd_grid(a):
    os.makedirs(a.out, exist_ok=True)
    pl = json.load(open(a.params))
    cont = model_moments(pl, n_workers=100_000, n_years_keep=36, n_age_cols=N_AGE, seed=987654)
    timer = ss_timer()
    rows = []
    zk, Qk = chain_from_files(os.path.join(KMV_CHAIN, "ygrid_combined.txt"),
                              os.path.join(KMV_CHAIN, "ymarkov_combined.txt"))
    t_kmv, A_kmv, err = timer(zk, Qk)
    rows.append({"chain": "KMV published 3x11 (reference)", "n_states": 33, "ss_seconds": t_kmv,
                 "ss_A": A_kmv, "ss_error": err})
    for n1, n2 in ((3, 11), (3, 15), (5, 11), (5, 15)):
        d = discretize_process(pl, n1=n1, n2=n2, match_var_log=True,
                               sim_kwargs=dict(n_workers=50_000, n_age_cols=N_AGE))
        m = moments_from_panel(simulate_discrete_panel(d, n_workers=100_000, n_years_keep=36),
                               n_age_cols=N_AGE)
        t, A, err = timer(d["combined"]["zgrid"], d["combined"]["Q"])
        r = {"chain": f"Poland v3 {n1}x{n2}", "n_states": n1 * n2, "ss_seconds": t, "ss_A": A,
             "ss_error": err, "ergodic_var_log_e": chain_props(d["combined"]["zgrid"], d["combined"]["Q"])["ergodic_var_log_e"]}
        for k in MOMENT_ORDER:
            r[f"chain_{k}"] = m[k]
            r[f"relgap_{k}"] = m[k] / cont[k] - 1
        rows.append(r)
        print(r, flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(a.out, "grid_size.csv"), index=False)
    json.dump({"continuous_25_55": cont, "rows": rows}, open(os.path.join(a.out, "grid_size.json"), "w"),
              indent=2, default=float)
    pd.set_option("display.width", 250)
    print(df[["chain", "n_states", "ss_seconds", "ergodic_var_log_e", "relgap_var_log_earns",
              "relgap_var_d1", "relgap_var_d5", "relgap_kurt_d1", "relgap_kurt_d5",
              "relgap_frac_d1_lt_10"]].to_string(index=False, float_format="%.4g"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["table", "grid"])
    ap.add_argument("--params", default="output/pol_v3/params.json")
    ap.add_argument("--chain", default="output/pol_v3")
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()
    {"table": cmd_table, "grid": cmd_grid}[a.what](a)


if __name__ == "__main__":
    main()
