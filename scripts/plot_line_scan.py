"""
Combines the original TikTak-to-DE line scan (a=0..1,
output/tiktak_line_scan.json) with its extension past the DE point (a=1..2,
output/tiktak_line_scan_extended.json) into one picture per country, with a
vertical line marking the DE (canonical) point at a=1 -- the plot form of
DECISIONS.md Sec 9's "the canonical point is the valley floor" result.
"""
from __future__ import annotations
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

orig = json.load(open("output/tiktak_line_scan.json"))
ext = json.load(open("output/tiktak_line_scan_extended.json"))

plt.rcParams.update({"figure.dpi": 150, "savefig.dpi": 150, "font.size": 10,
                     "axes.grid": True, "grid.alpha": 0.3})

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for ax, cc in zip(axes, ["FRA", "GER"]):
    a_orig = [p["a"] for p in orig[cc]["line_scan"]]
    f_orig = [p["f"] for p in orig[cc]["line_scan"]]
    a_ext = [p["a"] for p in ext[cc]["extension"]]
    f_ext = [p["f"] for p in ext[cc]["extension"]]

    ax.plot(a_orig, f_orig, "o-", color="#4472c4", lw=1.8, ms=5,
            label="TikTak -> DE (original scan)")
    ax.plot([1.0] + a_ext, [f_orig[-1]] + f_ext, "o-", color="#c0392b", lw=1.8, ms=5,
            label="DE -> extended (25-100% further)")
    ax.axvline(1.0, color="black", ls="--", lw=1.2, label="a=1: DE (canonical)")
    ax.set_yscale("log")
    ax.set_xlabel("a  (0 = TikTak point, 1 = DE/canonical, >1 = extended beyond DE)")
    ax.set_ylabel("objective f (log scale)")
    ax.set_title(f"{cc}: line scan, TikTak to DE and beyond")
    ax.legend(fontsize=8, loc="upper left")

fig.suptitle("The canonical (DE) estimate is the valley floor, not a point\n"
             "on a still-descending slope (DECISIONS.md Sec 9)")
fig.tight_layout(rect=(0, 0, 1, 0.92))
fig.savefig("output/pl_diagnostics/01_de_vs_tiktak_linescan.png")
print("Saved: output/pl_diagnostics/01_de_vs_tiktak_linescan.png")
