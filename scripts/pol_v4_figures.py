"""
Sec 32 (item 6): figures for Poland v4 (KMV specification), with KMV for
comparison. KMV side = KMV's paper parameters (E/earnings_estimation_paper,
km.KMV_PAPER_PARAMS) through the same simulator and design, and KMV's
published 3 x 11 chain (../annual-review/inputs/kmv_process).

Figures -> output/pol_v4/figures/:
  d1_distribution.png   1-year changes: CDF with GRID percentiles, density (log
                        scale), one-sided tails. Continuous processes in KMV's
                        design (one cross-section, 400k workers) and the chains
                        (ergodic start, quarterly expm(Q), annual flows summed)
  qq_vs_data.png        model vs GRID quantiles of 1-year changes
  impulse_response.png  E[z_j(t)] after a +1 sd jump (reset: decay beta + lambda),
                        and the autocorrelation of log e by lag for the
                        processes and the chains (= the expected response of log e
                        to a 1-unit deviation)
  var_changes.png       variance of 1- and 5-year changes: GRID by year with the
                        linear trend vs the model (v4 4-yr "d5" as fitted, v4
                        true 5-yr change, chains, v3)
  ergodic_chain.png     ergodic distribution of log e of the chains

Usage: python scripts/pol_v4_figures.py
"""

from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.ticker
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.linalg import expm
from scipy.stats import norm

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from kmv_earnings.simulate import simulate_log_annual_earnings, moments_from_panel
from kmv_earnings import kmv_method as km
import pol_v4 as v4
from pol_v3_figures import data_percentiles, PROBS

EST = "output/pol_v4/estimate"
OUT = "output/pol_v4/figures"
KMV_CHAIN = "../annual-review/inputs/kmv_process"
C_PL, C_KMV, C_PL2, C_V3 = "#c0392b", "#2c3e50", "#e67e22", "#8e44ad"


def proc_panel(p, n=400_000, years=6, seed=2468):
    s = dict(v4.BASE_SIM, n_workers=n, n_years_keep=years, seed=seed)
    s.pop("cross_section"); s.pop("d5_gap")
    return simulate_log_annual_earnings(p, **s)


def chain_panel(zgrid, Q, n=400_000, years=6, seed=97):
    """Ergodic start, quarterly P = expm(Q), annual earnings = sum of the 4
    quarterly levels exp(z), log."""
    P = expm(Q); P = np.clip(P, 0, None); P /= P.sum(1, keepdims=True)
    w, vv = np.linalg.eig(P.T)
    pi = np.abs(np.real(vv[:, np.argmin(np.abs(w - 1))])); pi /= pi.sum()
    rng = np.random.default_rng(seed)
    cum = np.cumsum(P, 1)
    s = np.searchsorted(np.cumsum(pi), rng.random(n))
    s = np.minimum(s, len(pi) - 1)
    out = np.empty((n, years)); flow = np.zeros(n)
    for q in range(4 * years):
        flow += np.exp(zgrid[s])
        if q % 4 == 3:
            out[:, q // 4] = np.log(flow); flow[:] = 0
        u = rng.random(n)
        s = np.minimum((u[:, None] > cum[s]).sum(1), len(pi) - 1)
    return out, P, pi


def chain_acf(zgrid, P, pi, lags_q):
    z = zgrid - pi @ zgrid
    v = pi @ z ** 2
    out, x = [], z.copy()
    Pk = np.eye(len(z))
    last = 0
    for L in lags_q:
        Pk = Pk @ np.linalg.matrix_power(P, L - last); last = L
        out.append(float((pi * z) @ (Pk @ z) / v))
    return np.array(out)


def main():
    os.makedirs(OUT, exist_ok=True)
    pl = json.load(open(os.path.join(EST, "params.json")))
    kmv = km.KMV_PAPER_PARAMS
    t = v4.load_targets()
    chains = {}
    for lab, d, c in (("Poland v4 chain 5x15", os.path.join(EST, "chain_ours_5x15"), C_PL2),
                      ("Poland v4 chain 3x11 (KMV procedure)", os.path.join(EST, "chain_kmv_3x11"), "#d35400"),
                      ("Poland v3 chain 5x15", "output/pol_v3", C_V3)):
        if os.path.exists(os.path.join(d, "income_process_Q.txt")):
            chains[lab] = (np.loadtxt(os.path.join(d, "income_process_zgrid.txt")),
                           np.loadtxt(os.path.join(d, "income_process_Q.txt")), c)
    chains["KMV published chain 3x11"] = (np.loadtxt(os.path.join(KMV_CHAIN, "ygrid_combined.txt")),
                                          np.loadtxt(os.path.join(KMV_CHAIN, "ymarkov_combined.txt")), C_KMV)
    ly = {"Poland v4 process": (proc_panel(pl), C_PL, "-"),
          "KMV process (paper params)": (proc_panel(kmv), "#7f8c8d", "-")}
    cpan = {}
    for lab, (zg, Q, c) in chains.items():
        y, P, pi = chain_panel(zg, Q)
        cpan[lab] = (y, P, pi, zg)
        if "v3" not in lab:
            ly[lab] = (y, c, "--")
    dq = data_percentiles()
    stats = {"data_quantiles": dq.tolist()}

    # ---- 1-year changes ----
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    for lab, (y, c, ls) in ly.items():
        d = np.sort(y[:, 1] - y[:, 0])
        F = (np.arange(len(d)) + 0.5) / len(d)
        ax[0].plot(d, F, color=c, ls=ls, lw=1.3, label=lab)
        h, e = np.histogram(d, bins=np.linspace(-4, 4, 161), density=True)
        ax[1].semilogy(0.5 * (e[1:] + e[:-1]), np.where(h > 0, h, np.nan), color=c, ls=ls, lw=1.1, label=lab)
        ax[2].semilogy(d, 1 - F, color=c, ls=ls, lw=1.3, label=f"{lab}: P(> x)")
        ax[2].semilogy(-d, F, color=c, ls=ls, lw=0.8, alpha=0.6)
        stats[lab] = {"var_d1": float(d.var()), "q": np.quantile(d, PROBS).tolist(),
                      "moments": moments_from_panel(y, cross_section=True, d5_gap=4)}
    ax[0].plot(dq, PROBS, "o", ms=5, mfc="white", mec="k", label="GRID Poland (mean pctl 2000-19)")
    ax[2].semilogy(dq[PROBS > 0.5], 1 - PROBS[PROBS > 0.5], "^", ms=5, mfc="white", mec="k",
                   label="GRID: P(change > x)")
    ax[2].semilogy(-dq[PROBS < 0.5], PROBS[PROBS < 0.5], "v", ms=5, mfc="white", mec="0.4",
                   label="GRID: P(change < -x) (thin lines: models)")
    x = np.linspace(-4, 4, 400)
    ax[1].semilogy(x, norm.pdf(x, scale=np.sqrt(t["var_d1"])), ":", color="0.5", label="normal, GRID variance")
    ax[0].set_xlim(-2, 2); ax[0].set_title("CDF, 1-year change"); ax[0].legend(fontsize=6.5)
    ax[1].set_ylim(1e-4, 20); ax[1].set_title("density (log scale)"); ax[1].legend(fontsize=6.5)
    ax[2].set_xlim(0, 4); ax[2].set_ylim(1e-5, 1); ax[2].set_title("one-sided tails (log scale)")
    ax[2].legend(fontsize=6.5)
    for a_ in ax:
        a_.set_xlabel("1-year change in log earnings")
    fig.suptitle("Poland v4 (KMV specification) vs KMV: 1-year changes, one ergodic cross-section", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "d1_distribution.png"), dpi=130); plt.close(fig)

    # ---- QQ ----
    fig, ax = plt.subplots(figsize=(5.6, 5.4))
    for lab, (y, c, ls) in ly.items():
        ax.plot(dq, stats[lab]["q"], "o", ls=ls, ms=3.5, lw=1, color=c, label=lab)
    lim = [min(dq.min(), -3), max(dq.max(), 3)]
    ax.plot(lim, lim, "-", color="0.75", lw=1, label="45 degrees")
    ax.set_xlabel("GRID Poland quantile (1-year change)"); ax.set_ylabel("model quantile")
    ax.set_title("model vs data quantiles, 1-year change (p1 ... p99.99)", fontsize=9)
    ax.legend(fontsize=6.5); fig.tight_layout(); fig.savefig(os.path.join(OUT, "qq_vs_data.png"), dpi=130)
    plt.close(fig)

    # ---- impulse responses + autocorrelation ----
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.4))
    yrs = np.linspace(0, 40, 801)
    for lab, p, c in (("Poland v4", pl, C_PL), ("KMV (paper params)", kmv, C_KMV)):
        for j, a_ in ((1, ax[0]), (2, ax[1])):
            b, s, lam = p[f"beta{j}"], p[f"sigma{j}"], p[f"lambda{j}"]
            a_.plot(yrs, s * np.exp(-(b + lam) * 4 * yrs), color=c, lw=1.6,
                    label=f"{lab}: sigma={s:.2f}, half-life {np.log(2)/(b+lam)/4:.2g} yr")
        d = v4.derived(p)
        acf = (d["var1"] * np.exp(-(p["beta1"] + p["lambda1"]) * 4 * yrs)
               + d["var2"] * np.exp(-(p["beta2"] + p["lambda2"]) * 4 * yrs)) / d["var_total"]
        ax[2].plot(yrs, acf, color=c, lw=1.6, label=f"{lab} process")
    lag_y = np.arange(0, 41)
    stats["acf_chains"] = {}
    for lab, (y, P, pi, zg) in cpan.items():
        a = chain_acf(zg, P, pi, list(4 * lag_y))
        stats["acf_chains"][lab] = a.tolist()
        ax[2].plot(lag_y, a, "o--", ms=2.5, lw=0.9, color=chains[lab][2], label=lab)
    ax[0].set_xlim(0, 2); ax[0].set_title("transitory: E[z1(t)] after a +1 sd jump (reset)")
    ax[1].set_title("persistent: E[z2(t)] after a +1 sd jump (reset)")
    ax[2].set_title("autocorrelation of quarterly log e at lag (years)")
    for a_ in ax:
        a_.set_xlabel("years after the shock"); a_.legend(fontsize=6.5); a_.axhline(0, color="0.7", lw=0.6)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "impulse_response.png"), dpi=130); plt.close(fig)

    # ---- variance of 1- and 5-year changes: model vs data ----
    panel = pd.read_csv("output/pol_v3/year_panel_pol_v3.csv")
    rep = json.load(open(os.path.join(EST, "report.json")))
    y5 = ly["Poland v4 process"][0]
    true5 = float(np.var(y5[:, 5] - y5[:, 0]))
    gap4 = float(np.var(y5[:, 4] - y5[:, 0]))
    v3 = json.load(open("output/pol_v3/report.json"))
    v3m = v3.get("oos_moments") or v3.get("canonical", {}).get("oos_moments", {})
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.4))
    for a_, col, lab in ((ax[0], "var_d1", "1-year"), (ax[1], "var_d5", "5-year")):
        s = panel[["year", col]].dropna()
        a_.plot(s.year, s[col], "o-", color="k", ms=4, lw=1, label=f"GRID Poland, var of {lab} change")
        b = np.polyfit(s.year, s[col], 1)
        a_.plot(s.year, np.polyval(b, s.year), ":", color="k", lw=1, label=f"linear trend (mean {s[col].mean():.3f})")
        if col == "var_d1":
            vals = [("v4 model", rep["oos_moments"]["var_d1"], C_PL, "-")]
        else:
            vals = [("v4 model, 4-yr gap (as fitted)", gap4, C_PL, "-"),
                    ("v4 model, true 5-yr change", true5, C_PL, "--")]
        for lab_c, (yc, *_r) in cpan.items():
            g = 1 if col == "var_d1" else 5
            vals.append((f"{lab_c}, {g}-yr", float(np.var(yc[:, g] - yc[:, 0])), chains[lab_c][2], "-."))
        if v3m:
            vals.append(("v3 model (lifecycle, ages 25-55)", v3m[col], C_V3, ":"))
        vals.append(("KMV target" if col == "var_d1" else "KMV target (4-yr gap)",
                     km.KMV_TARGETS[col], C_KMV, ":"))
        for lab_v, val, c, ls in vals:
            a_.axhline(val, color=c, ls=ls, lw=1.3, label=f"{lab_v}: {val:.3f}")
        stats[f"var_change_{col}"] = {l: v for l, v, *_ in vals}
        a_.set_title(f"variance of {lab} changes in log earnings"); a_.set_xlabel("base year")
        a_.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
        a_.legend(fontsize=6.5)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "var_changes.png"), dpi=130); plt.close(fig)

    # ---- ergodic distributions ----
    fig, ax = plt.subplots(figsize=(7, 4.4))
    bins = np.linspace(-6, 4, 51)
    for lab, (y, P, pi, zg) in cpan.items():
        e = np.exp(zg); e /= pi @ e; le = np.log(e)
        h, _ = np.histogram(le, bins=bins, weights=pi)
        v = float(pi @ (le - pi @ le) ** 2)
        ax.step(bins[:-1], h, where="post", color=chains[lab][2], lw=1.4, label=f"{lab}: var(log e) = {v:.3f}")
        stats[f"ergodic_var {lab}"] = v
    ax.set_xlabel("log e (mean e = 1), 0.2-wide bins"); ax.set_ylabel("ergodic mass")
    ax.set_title("ergodic distribution of the income chains"); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "ergodic_chain.png"), dpi=130); plt.close(fig)
    json.dump(stats, open(os.path.join(OUT, "figure_stats.json"), "w"), indent=2, default=float)
    print("saved figures to", OUT)


if __name__ == "__main__":
    main()
