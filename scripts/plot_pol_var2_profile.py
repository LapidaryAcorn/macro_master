"""
Plot OOS f against var2 for Poland's ergodic-variance profile (DECISIONS.md
Sec 23), including the canonical (free) estimate as the 5th point, and mark
the Monte-Carlo-noise band (+20% of the minimum f, per Sec 9's measured
~23% sd/mean at this objective scale) to show which tested var2 values are
statistically indistinguishable from the best one.
"""
from __future__ import annotations
import json
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

CANON_VAR2 = 0.5094416445471451
CANON_F = 0.00367  # Sec 18

summary = json.load(open("output/pol_var2_profile_summary.json"))
pts = [(CANON_VAR2, CANON_F, "canonical (free)")]
for v, r in summary["profiles"].items():
    pts.append((float(v), r["oos_f"], r["worst_moments"][0]["key"]))
pts.sort()

var2 = np.array([p[0] for p in pts])
f = np.array([p[1] for p in pts])
labels = [p[2] for p in pts]

f_min_profile = min(p[1] for p in pts if p[2] != "canonical (free)")
band_hi = f_min_profile * 1.20

fig, ax = plt.subplots(figsize=(8, 5.5))
ax.plot(var2, f, "o-", color="#1a9850", ms=8, lw=1.5, zorder=3)
for v, fv, lab in pts:
    ax.annotate(lab if lab == "canonical (free)" else lab,
               (v, fv), textcoords="offset points", xytext=(8, 6), fontsize=8)
ax.axhline(band_hi, color="grey", ls="--", lw=1,
          label=f"+20% of best pinned point ({f_min_profile:.4f}) -> {band_hi:.4f}\n"
                f"(Sec 9: Monte Carlo noise at this f-scale is ~23% of the mean)")
ax.axvline(CANON_VAR2, color="#1a9850", ls=":", lw=1, alpha=0.6)
ax.set_xlabel("var2 = lambda2 * sigma2^2 / (2 beta2)  (persistent-component ergodic variance)")
ax.set_ylabel("OOS objective f")
ax.set_yscale("log")
ax.set_title("Poland: OOS fit vs. the persistent component's ergodic variance\n"
            "(canonical/free estimate marked; grey line = noise-indistinguishable threshold)")
ax.legend(fontsize=8, loc="upper center")
fig.tight_layout()
os.makedirs("output/pl_diagnostics", exist_ok=True)
fig.savefig("output/pl_diagnostics/11_var2_profile.png", dpi=150)
print("Saved: output/pl_diagnostics/11_var2_profile.png")

print(f"\nf_min (best pinned point, var2=0.7) = {f_min_profile:.5f}")
print(f"+20%% noise-band threshold = {band_hi:.5f}")
print("Points within the noise band of the best pinned point:")
for v, fv, lab in pts:
    tag = "  <-- within band" if fv <= band_hi else ""
    print(f"  var2={v:.4f}  f={fv:.5f}  worst={lab}{tag}")
