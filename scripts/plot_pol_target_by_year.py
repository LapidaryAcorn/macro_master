"""
Sec 24: Poland's var_log_earns target by year -- v1 (GRID std_perm^2,
residualized log 3-year average) vs v2 (GRID std_loginc^2, annual log
earnings). The target is the simple mean over 2000-2020 of the yearly
variances (dashed lines). Writes output/pol_v2/target_by_year.png.
"""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

d = pd.read_csv("data/PL std_log_inc.csv")
d = d[(d.country == "POL") & (d.gender == "Male") & d.year.between(2000, 2020)]
v1 = json.load(open("targets/pol_grid.json"))["var_log_earns"]
v2 = json.load(open("targets/pol_grid_v2.json"))["var_log_earns"]

fig, ax = plt.subplots(figsize=(8.5, 5))
ax.plot(d.year, d.std_loginc ** 2, "o-", color="#1a9850", label="v2: std_loginc² (roczny log dochodu, surowy)")
ax.plot(d.year, d.std_perm ** 2, "s-", color="#7f7f7f", label="v1: std_perm² (log średniej 3-letniej, rezydualizowany)")
ax.axhline(v2, color="#1a9850", ls="--", lw=1)
ax.axhline(v1, color="#7f7f7f", ls="--", lw=1)
ax.text(2020.3, v2, f"cel v2 = {v2:.3f}", va="center", color="#1a9850", fontsize=9)
ax.text(2020.3, v1, f"cel v1 = {v1:.3f}", va="center", color="#7f7f7f", fontsize=9)
ax.annotate("", xy=(2019.5, v2), xytext=(2019.5, v1),
            arrowprops=dict(arrowstyle="<->", color="black", lw=0.8))
ax.text(2019.2, (v1 + v2) / 2, f"+{v2 - v1:.3f}", ha="right", va="center", fontsize=9)
ax.set_xlim(1999.5, 2022.8)
ax.set_xticks(range(2000, 2021, 2))
ax.set_ylabel("wariancja (std² z GRID)")
ax.set_title("Polska, mężczyźni 25–55: cel var_log_earns rok po roku\n"
             "(cel = zwykła średnia wariancji rocznych 2000–2020; std_perm brak dla 2000–2001)")
ax.legend(fontsize=8, loc="lower left")
fig.tight_layout()
fig.savefig("output/pol_v2/target_by_year.png", dpi=150)
print(d.assign(v_loginc=d.std_loginc**2, v_perm=d.std_perm**2)[["year", "v_loginc", "v_perm"]].round(3).to_string(index=False))
