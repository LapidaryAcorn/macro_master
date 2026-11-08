"""
Sec 28 (D): the Polish v3 process against the GRID data and against KMV
(Table 3 parameters, same simulator and design as ours: additive jumps,
lifecycle panel -- Sec 29.5 for why KMV's own published process differs).

Figures -> output/pol_v3/figures/:
  d1_distribution.png   1-year changes: model CDF with GRID percentiles, density
                        (log scale), tails (log survival of |change|)
  d5_distribution.png   5-year changes: same panels (GRID publishes no
                        percentiles of 5-year changes; the data enter as a normal
                        with GRID's variance, for scale only)
  qq_vs_data.png        model quantiles vs GRID quantiles of 1-year changes
  var_by_age.png        within-age variance of log annual earnings by age
  sample_paths.png      z1, z2 and z1 + z2 for a few workers (same random numbers)
  impulse_response.png  expected path of a component after a typical (1 sd) shock
  ergodic_chain.png     ergodic distribution of log e of the exported chains

Model moments are on ages 25-55 (31 panel columns), as in v3. Data percentiles
are the mean over 2000-2019 of GRID's annual percentiles (POL, men 25-55).

Usage: python scripts/pol_v3_figures.py [--params output/pol_v3/params.json]
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import norm

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from kmv_earnings.simulate import KMV_TABLE3_PARAMS, simulate_log_annual_earnings
from kmv_earnings.discretize import discretize_process

N_AGE = 31
OUT = "output/pol_v3/figures"
KMV_CHAIN = "../annual-review/inputs/kmv_process"
C_PL, C_KMV = "#c0392b", "#2c3e50"
PCOLS = ["p1", "p2_5", "p5", "p10", "p12_5", "p25", "p37_5", "p50", "p62_5", "p75",
         "p87_5", "p90", "p95", "p97_5", "p99", "p99_9", "p99_99"]
PROBS = np.array([1, 2.5, 5, 10, 12.5, 25, 37.5, 50, 62.5, 75, 87.5, 90, 95, 97.5, 99,
                  99.9, 99.99]) / 100


def data_percentiles():
    d = pd.read_csv("data/dane grid PL.csv")
    d = d[(d.year >= 2000) & (d.year <= 2019)]
    return np.array([d[f"{c}_resid_d1yr"].mean() for c in PCOLS])


def changes(logy, gap):
    y = logy[:, :N_AGE]
    return (y[:, gap:] - y[:, :-gap]).ravel()


def simulate_components(p, n=6, years=31, spq=6, seed=7):
    """Quarterly z1, z2 paths for n workers from z = 0 (same algorithm and RNG
    call order as simulate_log_annual_earnings, without time aggregation)."""
    rng = np.random.default_rng(seed)
    lam = np.array([p["lambda1"], p["lambda2"]]) / spq
    beta = np.array([p["beta1"], p["beta2"]]) / spq
    sig = np.array([p["sigma1"], p["sigma2"]])
    z = np.zeros((n, 2))
    out = np.empty((years * 4, n, 2))
    for s in range(years * 4 * spq):
        for j in range(2):
            u = rng.random(n); eps = rng.standard_normal(n) * sig[j]
            zd = z[:, j] * np.exp(-beta[j])
            z[:, j] = np.where(u < 1 - np.exp(-lam[j]), zd + eps, zd)
        if s % spq == spq - 1:
            out[s // spq] = z
    return out


def chain_logs(zgrid, P_or_Q, is_generator):
    from scipy.linalg import expm
    P = expm(P_or_Q) if is_generator else P_or_Q
    P = P / P.sum(1, keepdims=True)
    w, v = np.linalg.eig(P.T)
    pi = np.abs(np.real(v[:, np.argmin(np.abs(w - 1))])); pi /= pi.sum()
    e = np.exp(zgrid); e /= pi @ e
    return np.log(e), pi


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", default="output/pol_v3/params.json")
    ap.add_argument("--chain", default="output/pol_v3")
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--n-workers", type=int, default=100_000)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    pl = json.load(open(a.params))
    kmv = KMV_TABLE3_PARAMS
    sim = dict(n_workers=a.n_workers, n_years_keep=N_AGE, seed=2468)
    ly = {"Poland v3": simulate_log_annual_earnings(pl, **sim),
          "KMV Table 3": simulate_log_annual_earnings(kmv, **sim)}
    col = {"Poland v3": C_PL, "KMV Table 3": C_KMV}
    dq = data_percentiles()
    var_d5_data = json.load(open("targets/pol_grid_v3.json"))["var_d5"]
    stats = {}

    # ---- distributions of 1- and 5-year changes ----
    for gap, fname in ((1, "d1_distribution.png"), (5, "d5_distribution.png")):
        fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
        for lab, y in ly.items():
            d = np.sort(changes(y, gap))
            F = (np.arange(len(d)) + 0.5) / len(d)
            ax[0].plot(d, F, color=col[lab], lw=1.4, label=f"model: {lab}")
            h, e = np.histogram(d, bins=np.linspace(-4, 4, 161), density=True)
            ax[1].semilogy(0.5 * (e[1:] + e[:-1]), np.where(h > 0, h, np.nan), color=col[lab], lw=1.2, label=lab)
            ad = np.sort(np.abs(d)); S = 1 - (np.arange(len(ad)) + 0.5) / len(ad)
            ax[2].semilogy(ad, S, color=col[lab], lw=1.4, label=lab)
            stats[f"{lab}_d{gap}_var"] = float(d.var())
        if gap == 1:
            ax[0].plot(dq, PROBS, "o", ms=5, mfc="white", mec="k",
                       label="GRID Poland (mean percentiles 2000-19)")
            ax[2].semilogy(np.abs(dq[PROBS > 0.5]), 1 - PROBS[PROBS > 0.5], "^", ms=5,
                           mfc="white", mec="k", label="GRID, upper percentiles (1 - p)")
            ax[2].semilogy(np.abs(dq[PROBS < 0.5]), PROBS[PROBS < 0.5], "v", ms=5,
                           mfc="white", mec="0.4", label="GRID, lower percentiles (p)")
        x = np.linspace(-4, 4, 400)
        sd = np.sqrt(stats[f"Poland v3_d{gap}_var"] if gap == 1 else var_d5_data)
        ax[1].semilogy(x, norm.pdf(x, scale=sd), ":", color="0.5", label="normal, same variance")
        ax[0].set_xlim(-2, 2); ax[0].set_title(f"CDF, {gap}-year change"); ax[0].legend(fontsize=7)
        ax[1].set_ylim(1e-4, 20); ax[1].set_title("density (log scale)"); ax[1].legend(fontsize=7)
        ax[2].set_xlim(0, 4); ax[2].set_ylim(1e-5, 1)
        ax[2].set_title(f"tails: P(|change| > x) (log scale)"); ax[2].legend(fontsize=7)
        for x_ in ax:
            x_.set_xlabel(f"{gap}-year change in log earnings")
        fig.tight_layout(); fig.savefig(os.path.join(a.out, fname), dpi=130); plt.close(fig)

    # ---- QQ vs data ----
    fig, ax = plt.subplots(figsize=(5.6, 5.4))
    for lab, y in ly.items():
        mq = np.quantile(changes(y, 1), PROBS)
        ax.plot(dq, mq, "o-", ms=4, lw=1, color=col[lab], label=lab)
        stats[f"{lab}_qq_model_q"] = mq.tolist()
    lim = [min(dq.min(), -3), max(dq.max(), 3)]
    ax.plot(lim, lim, "--", color="0.6", lw=1, label="45 degrees")
    ax.set_xlabel("GRID Poland quantile (1-year change)"); ax.set_ylabel("model quantile")
    ax.set_title("model vs data quantiles, 1-year change\n(p1 ... p99.99)")
    ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(os.path.join(a.out, "qq_vs_data.png"), dpi=130)
    plt.close(fig)
    stats["data_quantiles"] = dq.tolist()

    # ---- variance by age (36-year panel, ages 25-60; shaded = outside GRID 25-55) ----
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for lab, p in (("Poland v3", pl), ("KMV Table 3", kmv)):
        y36 = simulate_log_annual_earnings(p, n_workers=a.n_workers, n_years_keep=36, seed=2468)
        v = y36.var(axis=0)
        ax.plot(25 + np.arange(36), v, color=col[lab], lw=1.6, label=f"{lab} (mean 25-55: {v[:N_AGE].mean():.3f})")
        stats[f"{lab}_var_by_age"] = v.tolist()
    ax.axvspan(55.5, 60.5, color="0.92", zorder=0)
    ax.axhline(json.load(open("targets/pol_grid_v3.json"))["var_log_earns"], color=C_PL, ls=":",
               lw=1, label="Poland target (std_resid², mean over ages)")
    ax.axhline(0.70, color=C_KMV, ls=":", lw=1, label="KMV target 0.70")
    ax.set_xlabel("age (model: entry at 25 with z = 0)"); ax.set_ylabel("var(log annual earnings)")
    ax.set_title("within-age variance of log earnings"); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "var_by_age.png"), dpi=130); plt.close(fig)

    # ---- sample paths ----
    fig, ax = plt.subplots(2, 3, figsize=(15, 6.5), sharex=True)
    t = 25 + np.arange(31 * 4) / 4
    for r, (lab, p) in enumerate((("Poland v3", pl), ("KMV Table 3", kmv))):
        z = simulate_components(p, n=5, seed=11)
        for w in range(5):
            ax[r, 0].plot(t, z[:, w, 0], lw=0.8)
            ax[r, 1].plot(t, z[:, w, 1], lw=1.0)
            ax[r, 2].plot(t, z[:, w, 0] + z[:, w, 1], lw=1.0)
        ax[r, 0].set_ylabel(lab)
    for c, tt in enumerate(("z1 (transitory)", "z2 (persistent)", "log earnings z1 + z2")):
        ax[0, c].set_title(tt); ax[1, c].set_xlabel("age")
    for c in range(3):
        lo = min(ax[0, c].get_ylim()[0], ax[1, c].get_ylim()[0]); hi = max(ax[0, c].get_ylim()[1], ax[1, c].get_ylim()[1])
        ax[0, c].set_ylim(lo, hi); ax[1, c].set_ylim(lo, hi)
    fig.suptitle("five workers, same random numbers in both rows (quarterly, from z = 0 at 25)")
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "sample_paths.png"), dpi=130); plt.close(fig)

    # ---- impulse responses ----
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
    yrs = np.linspace(0, 40, 801)
    for lab, p in (("Poland v3", pl), ("KMV Table 3", kmv)):
        for j, axx in ((1, ax[0]), (2, ax[1])):
            b, s, lam = p[f"beta{j}"], p[f"sigma{j}"], p[f"lambda{j}"]
            axx.plot(yrs, s * np.exp(-b * 4 * yrs), color=col[lab], lw=1.6,
                     label=f"{lab}: sigma={s:.2f}, half-life {np.log(2)/b/4:.2g} yr")
            if lab == "KMV Table 3":
                axx.plot(yrs, s * np.exp(-(b + lam) * 4 * yrs), "--", color=col[lab], lw=1,
                         label="KMV, reset-jump reading (Sec 29.5): decay beta + lambda")
    ax[0].set_xlim(0, 2); ax[0].set_title("transitory component: E[z1(t)] after a +1 sd jump")
    ax[1].set_title("persistent component: E[z2(t)] after a +1 sd jump")
    for x_ in ax:
        x_.set_xlabel("years after the shock"); x_.legend(fontsize=7); x_.axhline(0, color="0.7", lw=0.6)
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "impulse_response.png"), dpi=130); plt.close(fig)

    # ---- ergodic distribution of the chains ----
    fig, ax = plt.subplots(figsize=(7, 4.4))
    chains = {}
    if os.path.exists(os.path.join(a.chain, "income_process_Q.txt")):
        chains["Poland v3 chain (5x15)"] = (np.loadtxt(os.path.join(a.chain, "income_process_zgrid.txt")),
                                            np.loadtxt(os.path.join(a.chain, "income_process_Q.txt")), True, C_PL)
    d = discretize_process(kmv, match_var_log=True, sim_kwargs=dict(n_workers=50_000, n_age_cols=N_AGE))
    chains["KMV Table 3, our discretization (5x15)"] = (d["combined"]["zgrid"], d["combined"]["Q"], True, "#7f8c8d")
    chains["KMV published chain (3x11, annual-review)"] = (
        np.loadtxt(os.path.join(KMV_CHAIN, "ygrid_combined.txt")),
        np.loadtxt(os.path.join(KMV_CHAIN, "ymarkov_combined.txt")), True, C_KMV)
    bins = np.linspace(-6, 4, 51)
    for lab, (zg, M, gen, c) in chains.items():
        le, pi = chain_logs(zg, M, gen)
        h, _ = np.histogram(le, bins=bins, weights=pi)
        v = float(pi @ (le - pi @ le) ** 2)
        ax.step(bins[:-1], h, where="post", color=c, lw=1.5, label=f"{lab}: var(log e) = {v:.2f}")
        stats[f"ergodic_var_{lab}"] = v
    ax.set_xlabel("log e (mean e = 1), 0.2-wide bins"); ax.set_ylabel("ergodic mass")
    ax.set_title("ergodic distribution of the income chain"); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "ergodic_chain.png"), dpi=130); plt.close(fig)

    json.dump(stats, open(os.path.join(a.out, "figure_stats.json"), "w"), indent=2)
    print("saved figures to", a.out)


if __name__ == "__main__":
    main()
