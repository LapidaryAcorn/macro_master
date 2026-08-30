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

## S4. Swapping in our chain — the raw chains do not drop into an infinite-horizon model

**What was done.** `swap_chain.py <COUNTRY>` replaces the ARS income process with
`../kmv_grid_step1/output/<cc>/income_process_Q.txt` (the generator, so their
`expm(...)` call reconstructs our `P` — verified rows sum to 1, entries ≥ 0),
re-calibrates `(beta_hi, dbeta, omega)` to the **same** ARS targets (A = 20,
impact labor MPC = 0.20, SCF Lorenz curve) — the agreed "re-calibrate per chain"
convention — then rebuilds the GE steady state and reruns the monetary /
deficit IRFs and the Fig-3b decomposition. `results/swap_<cc>.json`.

**Headline result — the free-`β₂` chains carry too much ergodic earnings risk
for a stationary HANK.** The estimation matched GRID's *36-year working-life*
moments (§2). The ARS model is **infinite-horizon** and confronts households
with the chain's **ergodic** distribution, whose `var(log e)` is far higher
because a near-permanent component keeps accumulating dispersion the panel never
sees:

| chain | ergodic var(log e) | their KMV chain | data cross-section (GRID `std_log_inc²`) |
|---|---|---|---|
| GER (β₂ **pinned**, §8b) | 1.04 | 0.85 | 0.65 |
| FRA (β₂ free) | 1.04 | 0.85 | 0.48 |
| USA (β₂ free) | **1.79** | 0.85 | 0.96 |

**Calibration outcome:**

| | GER (pinned) | FRA (free) | USA (free) | baseline |
|---|---|---|---|---|
| A = 20 hit? | ✓ (20.00) | ✓ (20.00) | ✓ (20.05) | ✓ |
| SCF Lorenz hit? | ✓ (err 0.00) | ✓ (err −0.02) | ✗ (err −0.37) | ✓ |
| **MPC = 0.20 hit?** | **✓ (0.200)** | **✗ (0.028)** | **✗ (0.024)** | ✓ |
| `dbeta` (β-heterogeneity) | 0.070 (healthy) | 0.005 (collapsed) | ~0 (collapsed) | 0.088 |
| frac hand-to-mouth | 0.220 | 0.005 | 0.006 | 0.279 |
| β_ave shift vs baseline | +2% (0.954→0.975) | — (degenerate) | — (degenerate) | — |

For **USA and FRA the calibration degenerates**: the precautionary-saving motive
from the high ergodic risk means almost no households are hand-to-mouth, so the
aggregate MPC cannot be pushed to 0.20 by any `(β_hi, dβ, ω)` — the optimiser
collapses `dβ` to ~0 and lands at MPC ≈ 0.03. **This is structural, not an
optimiser artifact** — a 5-start (USA) / 5-start (FRA) multistart
(`_calib_multistart.py`, `results/calib_multistart.json`) converges to the same
degenerate point every time (USA cost 0.0848 ± 1e-6; FRA cost 0.0149). The
MPC = 0.20 target is genuinely unreachable with these chains in this model.
**This kills the KMV replication**:
with a low MPC the indirect (income) amplification is weak, and the monetary
decomposition shifts toward direct (year-1 direct share USA 0.39, FRA 0.46, vs
baseline 0.20). *(USA's decomposition is also unreliable — 5% residual from an
un-cleared asset market at the degenerate calibration.)*

**GER (β₂ pinned) is the one that works** — and its result is informative:

| GER swap | value | baseline | KMV Table 7 |
|---|---|---|---|
| β_ave | 0.975 | 0.954 | — |
| frac hand-to-mouth | 0.220 | 0.279 | — |
| MPC (labor / unweighted, impact) | 0.20 / 0.25 | 0.20 / 0.29 | — |
| monetary dY impact (HA) | 2.31 | 2.41 | — |
| monetary dY cum-40q (HA) | 24.2 | 28.3 | — |
| **decomposition (year 1): direct / indirect** | **30 / 70** | 20 / 80 | 19 / 80 |
| — indirect: labor / tax / cap-gains | 32 / 13 / 24 | 42 / 16 / 23 | 51 / 32 / −2 |

So even the well-behaved chain **shifts the monetary decomposition ~10 pp from
indirect toward direct** (80/20 → 70/30). "Indirect effects dominate" *survives*,
but is weaker than baseline and than KMV — mainly because our chain's
transitory/persistent structure passes less of the equilibrium wage response
through to consumption (indirect-labor share 42% → 32%).

### Diagnosis and the fix (supervisor decision)

Root cause is the **ergodic-vs-lifecycle-panel gap** (DECISIONS §2), realised in
the model. Two things need to change, both in the **chain export**, not the model:

1. **Drop `match_var_log` for the step-2 export.** It rescales the grid so the
   *36-year panel* var(log y) matches the continuous process — which, for a
   near-permanent component, *inflates* the ergodic variance further. An
   infinite-horizon model needs the chain's *ergodic* variance controlled.

2. **Rescale the combined grid so the chain's ergodic `var(log e)` hits a
   target.** Candidates:
   - **KMV's ≈ 0.85** — clean apples-to-apples with ARS ("same ergodic
     dispersion, our persistence/kurtosis shape"). Recommended for the
     *validation* exercise: it isolates what the process *shape* does.
   - **the data cross-section** (GRID `std_log_inc²`: USA 0.96, FRA 0.48,
     GER 0.65) — the GRID number *is* a stationary cross-section of prime-age
     males, so arguably the right target; but then the models differ in scale
     too and the calibration re-does more work.

3. Separately, **pinning `β₂` (§8b) helps USA** (ergodic var 1.79 → ~1.19) but
   **not FRA** (FRA's free `β₂ = 0.0052` is already *more* mean-reverting than
   the KMV anchor 0.0046 would give — pinning FRA moves the wrong way). FRA's
   problem is the export scale, not `β₂`. So the ergodic-var rescale (point 2)
   is the general fix; `β₂` pinning is a separate, USA/GER-only lever.

*A caveat on point 2.* FRA and GER have the **same** ergodic var(log e) (1.04)
yet FRA's calibration fails and GER's works. So the ergodic *variance* alone does
not determine whether a chain is usable — the process *shape* matters too. GER
has a higher persistent-shock arrival rate (`λ₂ = 0.0102` vs FRA 0.0070) and
more kurtosis, which puts more households near the borrowing constraint (frac
hand-to-mouth 0.22 vs FRA 0.005) and lets the MPC target be hit. The ergodic-var
rescale is necessary but may not be sufficient; the rebuilt-chain swap will show
whether it is.

**Recommendation:** rebuild the step-2 chains with an ergodic-`var(log e)`
target (≈ 0.85, matching ARS) instead of `match_var_log`; no re-estimation
needed. Then re-run the swaps. Do not proceed to the §8d var-Δ1 fix until the
scale question is settled — they interact.

> **Correction (see §S5).** Two things from the KMV paper revise this:
> (1) KMV's *own* implied ergodic variance at their Table 4 parameters is
> `λσ²/(2β) = 0.159 + 0.910 = 1.07`, **not** 0.85 — 0.85 is the ARS baseline
> *process*, a specific 33-state discretization, not KMV's target. So the
> rescale target, if any is used, is **1.07**, and it **only binds for USA**
> (1.79); FRA and GER (1.04) already match KMV almost exactly and need no
> rescaling. (2) KMV's death mechanism does *not* reduce cross-sectional
> earnings dispersion — newborns draw `z` from its ergodic distribution — so
> "the model uses the ergodic distribution" is right, but our FRA/GER chains
> are *already* at KMV's ergodic dispersion. The USA outlier is driven
> specifically by our `β₂ = 0.0066` (× `σ₂ = 1.32`) against KMV's `β₂ = 0.009`.

### Also settled here

- **`beta` moves little when the chain is well-scaled.** GER: β_ave +2%,
  β_lo +1%. So a well-behaved chain does *not* imply wildly different earnings
  risk — the apparent "big β move" for USA/FRA is the degenerate calibration,
  not a real signal.
- The **generator/`P` interface works** — `expm(income_process_Q.txt)`
  reconstructs our `P`, assertions pass, `n_e = 75` flows through.

---

## S5. Does KMV's stochastic death fix the one-asset calibration? No.

**Hypothesis (from the KMV paper).** KMV footnote 15: households die at
`ζ = 1/180` per quarter and newborns start at **zero wealth**; this is stated to
be the device that puts enough households at zero illiquid wealth to generate
the aggregate MPC. ARS's one-asset code is infinitely-lived with no such device
(**confirmed** — no death / mortality / rebirth / survival anywhere in
`annual-review/*.py` or the notebooks; SSJ has no built-in support either). So:
add it and see whether the USA/FRA calibrations stop degenerating.

**Implementation** (`death.py`): subclass SSJ's `HetBlock`, override
`make_endog_law_of_motion` with a death-modified policy lottery that moves a
fraction `ζ` of the mass at every `(β, e, a)` to asset index 0 each period;
discount the continuation by survival, `β → β(1−ζ)`; Blanchard–Yaari annuities
(`return → (1+r)/(1−ζ)`, so aggregate wealth is conserved) as the default, with
a no-annuity switch. Flows through both `steady_state` and `jacobian`. Newborns
keep their `(β, e)` state — but since the cross-sectional `(β, e)` distribution
is the stationary one at all times, newborns are ergodic draws, exactly as KMV.

**Result — death does not generate hand-to-mouth behaviour in the one-asset
model.** Homogeneous `β` + death + annuities, `β` calibrated to `A = 20`
(`_death_homog.py`, `results/death_homog.json`):

| chain | β | MPC (labor) | MPC (unweighted) | frac at `a = 0` |
|---|---|---|---|---|
| **baseline (their KMV chain)** | 0.993 | **0.026** | 0.043 | 0.008 |
| USA | 0.991 | 0.024 | 0.051 | 0.012 |
| FRA | 0.993 | 0.022 | 0.042 | 0.009 |
| GER | 0.993 | 0.021 | 0.041 | 0.012 |

MPC ≈ 0.02–0.03 for **every** chain, including the baseline KMV one (target
0.20). With ARS's `β`-heterogeneity *plus* death the calibration still
degenerates for USA (`dβ → 0`, MPC 0.02) and lands at MPC 0.05 for FRA with
`β_hi` pinned at 1.

**Why death works in KMV but not here.** KMV is **two-asset**. Death generates
MPC there because newborns start at zero *illiquid* wealth and rationally stay
there for years — adjusting illiquid holdings is costly — while holding liquid
wealth and running a high MPC ("wealthy hand-to-mouth"). In a **one-asset**
model there is no transaction cost pinning newborns at zero; a newborn with
average earnings saves out of poverty within a few quarters, so the stock at
zero wealth is ~1% and contributes almost nothing to the aggregate MPC.
**ARS's `β`-heterogeneity is precisely the one-asset substitute for KMV's
two-asset wealthy-hand-to-mouth mechanism** — a permanently-impatient type that
stays at the constraint. Bolting death onto the one-asset model adds an
ineffective mechanism, not a fix.

**Consequence for the diagnosis.** The USA/FRA calibration failure (§S4) is
*not* solved by the KMV death mechanism. The live options are:
1. **USA** — a targeted `β₂` fix. USA's ergodic var 1.79 vs KMV's 1.07 comes
   from `β₂ = 0.0066` combined with `σ₂ = 1.32`. Pin the ergodic variance of the
   persistent component to KMV's (DECISIONS §8b `ErgodicVarConstraint`) → `β₂ ≈
   0.0105`, ergodic var → ~1.2. Re-estimate USA (~50 min). FRA and GER need
   nothing on this axis (already ≈ KMV's 1.07).
2. **FRA — confirmed cannot calibrate in one-asset** (`_fra_final.py`,
   `results/fra_final.json`: 3 seeds incl. GER's working solution, all collapse
   to `dβ ≈ 0.005`, MPC 0.03). FRA's ergodic variance is fine (1.04 ≈ KMV's
   1.07); the failure is process *shape* — lower `λ₂ = 0.007` (vs GER 0.010) →
   persistent shocks rarer → households can self-insure them → few near the
   constraint. Probing wider `β`-heterogeneity: as `β_lo` falls, MPC rises to
   0.4–0.6 **but `A` collapses to 3–12 (never 20) and the SCF Lorenz goes badly
   off**. So (`A = 20`, MPC = 0.20, SCF Lorenz) are **jointly infeasible** for
   the FRA chain in the one-asset `β`-het family. GER has a `β_lo ≈ 0.92` sweet
   spot; FRA does not.

**Where this leaves step 2.** One-asset validates *as a mechanism*
(reproduces KMV's 20/80 decomposition, §S3) and works for GER. But:
  - **USA** needs a `β₂` re-estimation (pin ergodic var of the persistent
    component to KMV's 1.07 → `β₂ ≈ 0.0105`) — then it should calibrate.
  - **FRA** cannot be made to calibrate in one-asset by any lever short of
    changing the targets. The honest options: (a) report FRA with its
    model-implied MPC (~0.03–0.05) and weaker transmission, caveated; (b) use a
    French wealth distribution / MPC target instead of the US SCF ones; (c) go
    **two-asset**, where the liquid/illiquid transaction cost generates
    wealthy-hand-to-mouth regardless of the earnings-process shape.
  - **Relevance to Poland (the actual thesis target):** if the Polish earnings
    process turns out shaped like France's (low `λ₂`, rare persistent shocks),
    the one-asset model won't calibrate for Poland either. That is a genuine
    argument for building two-asset now rather than discovering the problem at
    the Polish stage. Supervisor decision.

**Files:** `death.py`, `_death_homog.py` / `results/death_homog.json`,
`_calib_death.py` / `results/calib_death.json` (ζ=1/180 no-annuity multistart,
also degenerate).

---

## Open (step 2)

- **The one-asset calibration with our chains (§S4–S5) — the fork.** GER works.
  USA needs a `β₂` re-estimation. **FRA cannot calibrate in one-asset at all**
  (confirmed). Decision: (a) one-asset, USA re-estimated, FRA reported caveated
  with a low MPC; (b) one-asset with country-specific wealth/MPC targets for
  FRA; or (c) build the **two-asset** model — which also de-risks Poland if the
  Polish process is France-shaped. **Blocks everything downstream.** Supervisor
  question. The KMV stochastic-death hypothesis was tested (§S5) and does not
  fix one-asset.
- ~~Chain export scale (rescale to 0.85)~~ — **withdrawn** (§S5 correction):
  0.85 is the ARS process, not KMV's target (1.07); FRA/GER already match, only
  USA is high and that is a `β₂` issue, not a rescale issue.
- **Whether the aggregate 80/20 decomposition match is enough**, or the
  channel-level composition (labor / tax / capital-gains split) also has to line
  up with KMV — the latter needs changes to the non-household blocks (asset
  pricing, fiscal rule). See §S3. Supervisor question.
- **Polish `beta` calibration target.** For USA/FRA/GER we calibrate to known
  asset targets (A = 20 = 500% of annual GDP, SCF Lorenz). Poland's wealth
  moments (NBP HFS / ECB HFCS) are not in hand, and Polish household wealth
  differs from the US in level *and* composition, so Polish `beta` will land
  elsewhere and that alone moves MPCs and transmission. "What do we calibrate
  Polish `beta` to, and what is the asset target" is a thesis-design question
  for Zoch, not an implementation detail.
- **S4 var Δ1 fix** (DECISIONS §8d) — after the scale question, discretize
  component 1 as a near-iid quarterly shock, re-check var Δ1 vs continuous and
  whether wealth/MPC/IRFs move.
