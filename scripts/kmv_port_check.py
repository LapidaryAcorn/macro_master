"""
Sec 29: does the simulation design of the KMV estimation code port
(Karibzhanov, Bank of Canada SDP 2020-13, github.com/ikarib/HANK_EarningsEstimation,
estimate.cu) reproduce KMV Table 2 at the Table 3 parameters?

Port design (read from estimate.cu): RESET jumps (z -> sigma*eps on arrival),
Euler decay (1 - delta/4) per step, 16 steps/year with arrival prob lambda/4,
start from N(0, sigma^2/(1 + 2 delta/lambda)) (the reset process's stationary
variance), 25-year burn-in, then 5 years; moments pooled over workers for one
year (level: year 0; changes: year 1 vs 0 and year 4 vs 0 -- sic, the port's
"5-year" change is log(y4/y0), a 4-year gap).

Compared with this repo's design (additive jumps, lifecycle panel from z = 0,
36 years, Sec 1-2) and with KMV's Table 2 model column.
"""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from kmv_earnings.simulate import KMV_TABLE3_PARAMS as P, model_moments, MOMENT_ORDER


def kurt(x):
    x = x - x.mean(); return float(np.mean(x**4) / np.mean(x**2) ** 2)


def port_moments(p, n=400_000, seed=3, five_gap=4):
    rng = np.random.default_rng(seed)
    lam = np.array([p["lambda1"], p["lambda2"]]); sig = np.array([p["sigma1"], p["sigma2"]])
    dl = np.array([p["beta1"], p["beta2"]])
    z = rng.standard_normal((n, 2)) * sig / np.sqrt(1 + 2 * dl / lam)
    nyear = max(5, five_gap + 1)
    zann = np.zeros((n, nyear))
    for t in range(-25, nyear):
        for q in range(16):
            r = rng.standard_normal((n, 2)); u = rng.random((n, 2))
            z = np.where(u > 1 - lam / 4, sig * r, (1 - dl / 4) * z)
            if t >= 0:
                zann[:, t] += np.exp(z.sum(1))
    ly = np.log(zann)
    d1 = ly[:, 1] - ly[:, 0]; d5 = ly[:, five_gap] - ly[:, 0]
    a = np.abs(d1)
    return {"var_log_earns": float(ly[:, 0].var()), "var_d1": float(d1.var()),
            "var_d5": float(d5.var()), "kurt_d1": kurt(d1), "kurt_d5": kurt(d5),
            "frac_d1_lt_10": float((a < .1).mean()), "frac_d1_lt_20": float((a < .2).mean()),
            "frac_d1_lt_50": float((a < .5).mean())}


def kmv_chain_jump_structure(path="../annual-review/inputs/kmv_process/ymarkov_combined.txt"):
    """KMV's own published 33-state generator: with ADDITIVE jumps the rate of
    jumping into state j depends on the origin state (normal centred on the
    origin); with RESET jumps it does not (lambda * w_j for every origin)."""
    Q = np.loadtxt(path)
    O = Q.copy(); np.fill_diagonal(O, 0)
    distinct = np.unique(np.round(O[O > 0], 5))
    per_col = [sorted(set(np.round(O[:, j][O[:, j] > 0], 5))) for j in range(len(Q))]
    return {"n_states": len(Q), "n_distinct_offdiag_rates": int(len(distinct)),
            "distinct_rates": distinct.tolist(),
            "max_distinct_rates_into_one_state": max(len(c) for c in per_col)}


if __name__ == "__main__":
    kmv_model = dict(zip(MOMENT_ORDER, [0.70, 0.23, 0.46, 16.5, 12.1, 0.56, 0.67, 0.85]))
    out = {"kmv_table2_model": kmv_model,
           "repo_additive_lifecycle": model_moments(P, n_workers=100_000),
           "port_reset_ergodic_gap4": port_moments(P, five_gap=4),
           "port_reset_ergodic_gap5": port_moments(P, five_gap=5)}
    print(f"{'moment':15s}" + "".join(f"{k[:22]:>24s}" for k in out))
    for m in MOMENT_ORDER:
        print(f"{m:15s}" + "".join(f"{out[k][m]:24.3f}" for k in out))
    out["kmv_chain_jump_structure"] = kmv_chain_jump_structure()
    print(out["kmv_chain_jump_structure"])
    json.dump(out, open("output/pol_v3/kmv_method/port_check.json", "w"), indent=2)
