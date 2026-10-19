"""
Sec 24: Poland's ergodic-variance profile on the v2 targets (var_log_earns
from GRID std_loginc), same layout as plot_pol_var2_profile.py (Sec 23):
OOS f vs var2, +20%-of-minimum noise band, Germany's curve overlaid.
The v1 Polish curve (std_perm target) is drawn faded for comparison.

Reads output/pol_v2/close_step1_report.json (canonical v2 free estimate) and
output/pol_v2/var2_profile_summary.json (pinned points, incl. the control
point at the canonical v2 var2). Writes output/pol_v2/11_var2_profile_v2.png;
the v1 figure (output/pl_diagnostics/11_var2_profile.png) is left untouched.
"""
from __future__ import annotations
import json
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

V1_CANON_VAR2, V1_CANON_F = 0.5094416445471451, 0.00367
GER_FREE_F = 0.0025192182680727334   # same OOS panel, see plot_pol_var2_profile.py

rep = json.load(open("output/pol_v2/close_step1_report.json"))
pc = rep["params_canonical"]
canon_var2 = pc["lambda2"] * pc["sigma2"] ** 2 / (2 * pc["beta2"])
canon_f = rep["oos_f"]["de" if rep["winner"] == "DE" else "tiktak"]

summary = json.load(open("output/pol_v2/var2_profile_summary.json"))
pts = [(canon_var2, canon_f, "canonical v2 (free)")]
for v, r in summary["profiles"].items():
    is_control = abs(float(v) - canon_var2) < 5e-4
    lab = "control (pinned at canonical)" if is_control else r["worst_moments"][0]["key"]
    pts.append((float(v), r["oos_f"], lab))
pts.sort()
var2 = np.array([p[0] for p in pts])
f = np.array([p[1] for p in pts])
band_hi = canon_f * 1.20

fig, ax = plt.subplots(figsize=(8.5, 5.5))
ax.plot(var2, f, "o-", color="#1a9850", ms=8, lw=1.5, zorder=3,
        label="Poland v2 (var_log from std_loginc)")
_used = []
for v, fv, lab in pts:
    dy = 6
    while any(abs(v - ov) < 0.03 and abs(dy - ody) < 10 for ov, ody in _used):
        dy -= 14
    _used.append((v, dy))
    ax.annotate(lab, (v, fv), textcoords="offset points", xytext=(8, dy), fontsize=8)
ax.axhline(band_hi, color="grey", ls="--", lw=1,
           label=f"+20% of the v2 minimum ({canon_f:.4f}) -> {band_hi:.4f}")
ax.axvline(canon_var2, color="#1a9850", ls=":", lw=1, alpha=0.6)

v1_path = "output/pol_var2_profile_summary.json"
if os.path.exists(v1_path):
    v1 = json.load(open(v1_path))
    v1_pts = sorted([(V1_CANON_VAR2, V1_CANON_F)]
                    + [(float(v), r["oos_f"]) for v, r in v1["profiles"].items()])
    ax.plot([p[0] for p in v1_pts], [p[1] for p in v1_pts], "o:", color="#1a9850",
            alpha=0.35, ms=5, lw=1, zorder=1, label="Poland v1 (std_perm target, Sec 23)")

ger_path = "output/ger_var2_profile_summary.json"
if os.path.exists(ger_path):
    g = json.load(open(ger_path))
    gp = sorted([(g["free_var2"], GER_FREE_F)]
                + [(float(v), r["oos_f"]) for v, r in g["profiles"].items()])
    ax.plot([p[0] for p in gp], [p[1] for p in gp], "s--", color="#c0392b",
            ms=7, lw=1.5, zorder=2, label="Germany")

ax.set_xlabel("var2 = lambda2 * sigma2^2 / (2 beta2)  (persistent-component ergodic variance)")
ax.set_ylabel("OOS objective f")
ax.set_yscale("log")
ax.set_title("Poland (v2 targets) vs Germany: OOS fit vs. the persistent "
             "component's ergodic variance\n(grey line = noise-indistinguishable "
             "threshold; faded = Poland v1)")
ax.legend(fontsize=8, loc="best")
fig.tight_layout()
fig.savefig("output/pol_v2/11_var2_profile_v2.png", dpi=150)
print("Saved: output/pol_v2/11_var2_profile_v2.png")
print(f"canonical v2 var2={canon_var2:.4f}  f={canon_f:.5f}  band={band_hi:.5f}")
for v, fv, lab in pts:
    print(f"  var2={v:.4f}  f={fv:.5f}  {lab}{'  <-- within band' if fv <= band_hi else ''}")
