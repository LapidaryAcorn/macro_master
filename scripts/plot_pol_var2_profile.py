"""
Plot OOS f against var2 for Poland's ergodic-variance profile (DECISIONS.md
Sec 23), including the canonical (free) estimate, and mark the Monte-Carlo-
noise band (+20% of the MINIMUM f -- the canonical free estimate, not the
best PINNED point, corrected 2026-10-12) to show which tested var2 values
are statistically indistinguishable from the actual optimum.

Optionally overlays:
  - the var2=0.5094 control point (same pipeline as the profile, but pinned
    at exactly the canonical's own var2 -- should reproduce the canonical
    free estimate if the pin is genuinely non-binding there), once
    output/pol_var2_0.5094/profile_report.json exists;
  - Germany's own profile (output/ger_var2_profile_summary.json), once it
    exists, for a direct identification-strength comparison on one plot.
"""
from __future__ import annotations
import json
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

CANON_VAR2 = 0.5094416445471451
CANON_F = 0.00367  # Sec 18, Poland's canonical free-estimate OOS f -- the true minimum

summary = json.load(open("output/pol_var2_profile_summary.json"))
pts = [(CANON_VAR2, CANON_F, "canonical (free)")]
for v, r in summary["profiles"].items():
    pts.append((float(v), r["oos_f"], r["worst_moments"][0]["key"]))

control_path = "output/pol_var2_0.5094/profile_report.json"
have_control = os.path.exists(control_path)
if have_control:
    c = json.load(open(control_path))
    pts.append((0.5094, c["oos_f"], "control (pinned 0.5094)"))

pts.sort()
var2 = np.array([p[0] for p in pts])
f = np.array([p[1] for p in pts])

# Noise band: +20% of the MINIMUM (the canonical free estimate), not the
# best PINNED point -- the pinned grid are all constrained re-estimates of a
# DIFFERENT (5-parameter) problem, so the true reference optimum is the free
# 6-parameter estimate itself.
f_min = CANON_F
band_hi = f_min * 1.20

fig, ax = plt.subplots(figsize=(8.5, 5.5))
ax.plot(var2, f, "o-", color="#1a9850", ms=8, lw=1.5, zorder=3, label="Poland")
for v, fv, lab in pts:
    ax.annotate(lab, (v, fv), textcoords="offset points", xytext=(8, 6), fontsize=8)
ax.axhline(band_hi, color="grey", ls="--", lw=1,
          label=f"+20% of the minimum ({f_min:.4f}) -> {band_hi:.4f}\n"
                f"(Sec 9: Monte Carlo noise at this f-scale is ~23% of the mean)")
ax.axvline(CANON_VAR2, color="#1a9850", ls=":", lw=1, alpha=0.6)

ger_path = "output/ger_var2_profile_summary.json"
if os.path.exists(ger_path):
    ger_summary = json.load(open(ger_path))
    ger_free_var2 = ger_summary["free_var2"]
    ger_pts = [(ger_free_var2, 0.00220, "GER canonical (free)")]  # Sec 21
    for v, r in ger_summary["profiles"].items():
        ger_pts.append((float(v), r["oos_f"], None))
    ger_pts.sort()
    gv = [p[0] for p in ger_pts]
    gf = [p[1] for p in ger_pts]
    ax.plot(gv, gf, "s--", color="#c0392b", ms=7, lw=1.5, zorder=2, label="Germany")

ax.set_xlabel("var2 = lambda2 * sigma2^2 / (2 beta2)  (persistent-component ergodic variance)")
ax.set_ylabel("OOS objective f")
ax.set_yscale("log")
title = "OOS fit vs. the persistent component's ergodic variance"
if os.path.exists(ger_path):
    title = "Poland vs Germany: " + title
else:
    title = "Poland: " + title
ax.set_title(title + "\n(canonical/free estimate marked; grey line = noise-indistinguishable threshold)")
ax.legend(fontsize=8, loc="best")
fig.tight_layout()
os.makedirs("output/pl_diagnostics", exist_ok=True)
fig.savefig("output/pl_diagnostics/11_var2_profile.png", dpi=150)
print("Saved: output/pl_diagnostics/11_var2_profile.png")

print(f"\nf_min (canonical free estimate) = {f_min:.5f}")
print(f"+20% noise-band threshold = {band_hi:.5f}")
print("Points within the noise band of the TRUE minimum:")
for v, fv, lab in pts:
    tag = "  <-- within band" if fv <= band_hi else ""
    print(f"  var2={v:.4f}  f={fv:.5f}  worst={lab}{tag}")
if not have_control:
    print(f"\n(control point not yet computed -- run "
          f"`python scripts/pol_beta2_profile.py --var2 0.5094` and re-run this script)")
if not os.path.exists(ger_path):
    print("(GER profile not yet computed -- run scripts/ger_beta2_profile.py and re-run this script)")
