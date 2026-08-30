# Step 2 log — plugging the estimated chain into a one-asset HANK

Companion to `../kmv_grid_step1/DECISIONS.md` (see its §16 for the interface
notes this builds on). Same format: **what → why → how verified → implication.**

Scope reminder (DECISIONS.md §0): the thesis target is **Poland**. USA/FRA/GER
here are validation — evidence the pipeline reproduces KMV's mechanism on US
data before it is applied to Poland.

The step-2 model is **Auclert–Rognlie–Straub (2025, *Annual Review of
Economics*)**, `shade-econ/annual-review` — a compact one-asset HANK with
discount-factor heterogeneity, built on the sequence-space Jacobian toolkit.
Cloned to `../annual-review`.

---

## S1. Baseline reproduced unchanged

**What.** `baseline.py` reproduces the ARS one-asset HANK — their own KMV income
process (`inputs/kmv_process/*`, 33 states), their calibration (`hh_params.json`),
their GE blocks (copied verbatim from `Annual Review main.ipynb`) — and records
every number in `results/baseline.json`. Everything later is judged against this.

**Environment.** Python 3.14.6 (only interpreter on the machine). PyPI
`sequence-jacobian` is 1.0.0 (2022) and too old for the notebooks; installed from
**GitHub master** instead — imports and runs fine, `numba 0.67.0` jits on 3.14.
Also installed `nbconvert`, `jupytext`, `pypdf`, `ipykernel`.

**How verified.** Executed their `Annual Review main.ipynb` end-to-end headless
(`nbconvert --execute`, no errors) and diffed every key quantity against
`baseline.py`: steady state, monetary IRFs (HA/TA/RA), the Fig 3b decomposition,
the deficit-tax-cut IRFs. **Identical to full floating-point precision** (TA/RA
deficit impact differ at 1e-15 — FP noise). So the GE-block transcription is
exact and the baseline is faithful.

**Baseline numbers (their KMV chain, their calibration).**

| | value | check |
|---|---|---|
| steady state A / C (all models) | 20.0 / 0.80 | asset & goods mkt, NKPC residual ≈ 0 |
| HA frac at borrowing constraint | 0.279 | ~28% hand-to-mouth |
| HA wealth: top-10% / top-1% share | 0.753 / 0.223 | SCF 2019 ≈ 0.76 / 0.34 — top tail undershot (known one-asset limitation) |
| aggregate MPC (labor, impact) | 0.200 | calibration target, hit exactly |
| aggregate MPC (unweighted, impact) | 0.289 | matches their "≈ 0.29" |
| share of transfer spent in year 1 | 0.470 | matches their "slightly below 50%" |
| monetary shock, impact dY: HA / TA / RA | 2.41 / 2.21 / 1.99 | HA amplifies RA by ~21% on impact |
| monetary shock, cumulative-40q dY: HA / TA / RA | 28.3 / 21.8 / 19.6 | HA ~44% above RA cumulative |
| deficit tax cut, impact dY: HA / TA / RA | 0.52 / 0.21 / 0.00 | RA ≈ Ricardian, HA large |

`results/baseline_irfs.png` plots the three panels.

---

## S2. Frequency — **quarterly**, settled

**What.** The ARS model is quarterly. Our exported chain is quarterly
(`income_process_P.txt` = one-quarter transition), so it is fed **directly, no
`P⁴`**.

**How verified (multiple independent signals, all consistent).**
- Calibration notebook: `A_target = 20 # target 500% of annual GDP assets, or
  2000% of quarterly GDP`; `mpc_target = 0.20 # target quarterly MPC`.
- Every IRF plot: `plt.xlabel('Quarter')`; `T = 26` / `T = 400` quarters.
- "average **quarterly** MPC unweighted ≈ 0.29"; "slightly below 50% spent in
  the **first year**" = sum of the first 4 periods.
- `kappa = 0.01 # quarterly NKPC slope`.
- `household.py`: "asset income on A = 20 at r = 0.005 is 0.1" ⇒ with Y quarterly,
  A / annual-GDP = 5 ✓ (US ratio). `r = 0.005`/qtr = ~2%/yr.
- Their `ymarkov_combined.txt` is the KMV **quarterly** generator.

**Implication.** No frequency adjustment. The §8d fix (component 1 decays faster
than a quarterly chain resolves) stays relevant and will be applied at step S4.

---

## S3. One-asset vs two-asset — can the one-asset model reproduce KMV's headline?

**Question (raised, not to decide alone).** KMV's decomposition rests on the
*wealthy hand-to-mouth* structure that two liquid/illiquid assets generate. We
start one-asset. Does the one-asset model reproduce the qualitative result at
all, or does replicating the decomposition require two assets?

**Finding — the one-asset model reproduces the headline split almost exactly.**
The ARS one-asset model (β-heterogeneity instead of two assets) gives, for the
monetary shock, first-year-average consumption decomposition:

| component | ARS one-asset | KMV (2018) Table 7 col 1 |
|---|---|---|
| **direct** (interest rate) | **19.7%** | **19%** |
| **indirect total** | **80.3%** | **80%** |
| — indirect: labor income (w) | 41.8% | 51% |
| — indirect: transfers / taxes (T) | 16.0% | 32% |
| — indirect: asset returns / capital gains | 22.5% | −2% |

The **20 / 80 direct–indirect split is reproduced**. The composition of the
*indirect* effects differs, and the differences trace to the **non-household
blocks, not to one- vs two-asset**:
- ARS's `capitalization` block lets the share price rise with the rate cut →
  large positive capital-gains channel. KMV deliberately neutralise
  countercyclical profits (φ = ω̄) → that channel ≈ 0.
- ARS's fiscal rule ("constant real amount owed next period", a fast fiscal
  response) → smaller transfer channel than KMV's "lump-sum transfers adjust".

**Assessment.** For the validation milestone — reproduce KMV's finding that
*indirect GE effects dominate monetary transmission* — **one-asset with a
matched aggregate MPC is sufficient**, and does it quantitatively on the 20/80
headline. Two assets would be needed for: the cross-sectional MPC / consumption-
response distribution (KMV Fig 5–6), the wealthy-HtM mechanism *as such*, the
portfolio-rebalancing sub-channel of the direct effect, and the very top wealth
tail (one-asset undershoots top-1% share, 0.22 vs 0.34). None of those is the
headline decomposition. **Recommendation: proceed on the one-asset path; revisit
two-asset only if the thesis needs the distributional cross-section or the
portfolio channel.** Supervisor's call.

---

## Open (step 2)

- **`beta` calibration convention.** `beta` (the discount-factor grid) is
  calibrated to an asset target (A = 20) *given the income process*. Swapping our
  chain changes the calibrated `beta`. Any baseline-vs-ours comparison needs a
  stated convention: (a) re-calibrate `beta` per chain (matches wealth by
  construction, isolates the effect of the chain on *dynamics* / MPCs / IRFs), or
  (b) hold `beta` fixed at the baseline value (lets the chain move the wealth
  distribution too). (a) is the cleaner validation design. Supervisor question.
- **S4 var Δ1 fix** — discretize component 1 as a near-iid quarterly shock
  (DECISIONS §8d), then re-check chain var Δ1 vs continuous and whether the
  wealth/MPC picture moves vs S3.
