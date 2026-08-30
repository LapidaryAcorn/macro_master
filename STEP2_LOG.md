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

**But death + `β`-heterogeneity *together* DO rescue FRA** (`_death_betahet.py`,
`results/death_betahet.json`). Re-calibrating *both* mechanisms jointly to the
three targets:

| chain | converges? | MPC | A | `dβ` | `β_lo` | frac at 0 |
|---|---|---|---|---|---|---|
| FRA | **yes, exactly** (resid ~1e-12) | 0.20 | 20.0 | 0.093 | 0.905 | 0.225 |
| baseline | no — overshoots, `dβ` → 0.40 bound, MPC 0.35 | — | 20.0 | 0.40 | 0.60 | 0.354 |
| GER | (running) | | | | | |

So the mechanisms are **complementary, not redundant**: `β`-het supplies the
permanently-impatient type, death supplies the zero-wealth injection, and for
FRA's low-`λ₂` chain **you need both** to get enough households to the
constraint — `β`-het alone (§S4) left FRA stuck at MPC 0.03. The *baseline*
chain already has enough hand-to-mouth from `β`-het alone, so adding death there
overshoots (that solve hit the `dβ` bound and is not a clean calibration — a
lower `dβ` should work; the multistart bounds were too wide).

**Revised reading:** death is *not* inert. The right one-asset structure is
**ARS's `β`-heterogeneity plus KMV's stochastic death** — what KMV actually have
merged with what ARS actually have. Full swap re-run (`swap_death.py`,
`results/swap_death_*.json`; death block + joint `(β_hi, dβ, ω)` calibration,
then steady state / wealth / MPC / IRFs / decomposition vs the infinitely-lived
baseline):

| with death + β-het (Jacobian-fixed) | calibrates? | β_ave | frac HtM | MPC | mon. dY impact (cum-40) | decomp yr1 direct / indirect |
|---|---|---|---|---|---|---|
| baseline (their KMV chain) | — | 0.954 | 0.279 | 0.20 | 2.41 (28.3) | 20 / 80 |
| **FRA** | **yes, exact** | 0.970 (+2%) | 0.225 | 0.200 | 2.29 (−5%) / 23.6 (−17%) | **25.5 / 74.5** |
| **GER** | **yes, exact** | 0.976 (+2%) | 0.250 | 0.200 | 2.33 (−4%) / 24.0 (−15%) | **24 / 76** |
| **USA** (β₂-pinned re-estimate) | **no** | 0.990 | 0.011 | 0.029 | 1.85 (−23%) | degenerate |

*(indirect breakdown, FRA / GER: labour 33 / 34.5, tax 14 / 14.5, capital-gains
27.5 / 27 — of the total, in %.)*

**FRA and GER work in one-asset** — clean calibration (residuals ~1e-12), MPC
0.20, monetary IRF impact within 5% of baseline, **"indirect dominates"
reproduced: ≈75% indirect, vs KMV's 81% and baseline's 80%.** The direct share
(~25%) runs a few points above KMV's 19%, because our estimated transitory/
persistent structure passes less of the equilibrium wage response through to
consumption (indirect-labour share 33% vs baseline's 42%). β moves only +2%.
This **makes fork (c) two-asset unnecessary for FRA/GER**.

> **Jacobian bug fixed** (`death.py`, 2026-08-31). The `ForwardShockableDeathLottery1D`
> transition inherited the plain `PolicyLottery1D` `forward` / `expectation`
> instead of the death-modified ones — `het_block` builds `law_of_motion` from
> the *shockable* transitions and propagates the fake-news curlyE / curlyD
> vectors through them, so every Jacobian row at horizon ≥ 1 ignored death (the
> impact row stayed right). Symptom: the decomposition identity was ~10% off and
> `Jacobian` vs `impulse_nonlinear` diverged 45% at later horizons. After the fix
> (multiple-inherit from `DeathLottery1D`), the identity holds to 1e-6 and
> Jac-vs-nonlinear is 0.8% (linearization error). The table above and all
> `swap_death_*.json` are post-fix. `_decomp_jac_check.py` is the regression test.

**USA still does not calibrate, even with β₂ pinned *and* death.** The
β₂-pin re-estimate (DECISIONS §8b, now applied to USA — `β₂ = 0.0136`, ergodic
sd of the persistent component = KMV's 0.954) brings the *analytic* ergodic var
to ~1.05, but the exported chain lands at 1.19 (the `match_var_log` rescale
re-inflates it), and rescaling it down to 1.04 (= FRA/GER) or 1.07 (= KMV) still
does not fix it — `dβ → 0`, MPC 0.03. **The block is process *shape*, not
variance:** USA's GRID chain has `kurt Δ1y = 12.8` against KMV's 17.8 (GRID-USA
is a more volatile *and less leptokurtic* sample, DECISIONS §12). Without the
leptokurtic tail — the rare large shocks — too few households ever reach the
constraint, whatever the discount factor.

**Where this leaves the fork (much narrower now):**
- **FRA, GER:** one-asset + death. Done, works. No two-asset needed.
- **USA:** the clean answer is that **the S1 baseline (KMV's own chain) *is* the
  US validation** — it reproduces KMV's 20/80 decomposition (§S3). Our
  GRID-USA chain is then a *sensitivity* ("what if US earnings were as GRID
  measures them, less leptokurtic and more volatile"), and its answer — the MPC
  and transmission collapse — is itself a result about how much the tail of the
  earnings process matters. It does not need to "calibrate".
- **Poland:** still the open risk. If the Polish process is USA-shaped (thin
  tails), one-asset + death won't calibrate for Poland either. If it is
  FRA/GER-shaped, it will. Two-asset remains the hedge — see §S6 cost (~2 wk).

**Decomposition identity — fixed** (was ~10% off with the death block; now
1e-6). See the Jacobian-bug note above. The 25/75 figures are final, not
provisional.

**Consequence for the diagnosis (pre-`_death_betahet` — superseded above for
FRA).** The USA/FRA calibration failure (§S4) is *not* solved by death *alone*.
The live options were:
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
### The Poland argument for building two-asset now

**The thesis target is Poland.** The Polish earnings process is not yet estimated
(data pending). If it comes out **France-shaped** — low `λ₂`, rare persistent
shocks — then, by exactly the mechanism confirmed for FRA above, the **one-asset
model will not calibrate for Poland**: `(A, MPC, wealth distribution)` will be
jointly infeasible, and we would discover this only *after* building the entire
monetary/fiscal analysis on the one-asset path. France is not an unusual case —
lower `λ₂` than the US is a plausible European pattern, and Poland's labour
market (large informal sector, high job turnover at the low end, compressed
formal wage ladder) could go either way. Committing to one-asset is a bet that
Poland is Germany-shaped, taken before the Polish data are in.

**KMV's own argument reinforces this.** KMV (2018) build two assets precisely
because a **liquid-wealth-only** calibration "abstracts entirely from capital",
and the quantity and price responses of capital are a material part of the
**indirect** channel — which is the exact decomposition this thesis replicates.
A one-asset model with equity + bonds (ARS's structure, §S1) has a thin version
of this (the `capitalization` block), but no investment margin and no capital
stock; KMV's Table 2 indirect effects run partly through capital that a
one-asset model does not have. So even where one-asset *calibrates* (GER), the
channel decomposition it produces is structurally incomplete relative to what
the thesis is trying to reproduce (this is the §S3 "composition differs"
finding, seen from the model-structure side).

**Files:** `death.py`, `_death_homog.py`/`results/death_homog.json`,
`_death_betahet.py`/`results/death_betahet.json` (death + β-het, both
re-calibrated), `_calib_death.py`, `_fra_betahet.py`, `_fra_final.py`.

---

## S6. Costing the two-asset option (fork choice (c))

**Does a two-asset base exist? Yes, and it runs.** `sequence_jacobian.examples.
two_asset` is a complete, tested two-asset HANK (Auclert–Bardóczy–Rognlie–Straub
2021, Econometrica — the Kaplan–Violante/KMV structure): liquid `b`, illiquid
`a`, convex adjustment cost `Ψ(a', a; χ₀, χ₁, χ₂)`, capital + investment,
sticky prices *and* wages. `two_asset.dag()` solves its steady state in **6 s**
on this machine (`β = 0.973`, `χ₁ = 4.88`, illiquid `A = 13`, liquid `B = 1.04`,
asset market clears). The hard part — the 3-D `(z, b, a)` EGM het block with the
adjustment cost — is **already written** (`hetblocks/hh_twoasset.py`).

**The chain swap is identical to what we have done.** The z-process enters
through `make_grids` → `markov_rouwenhorst(rho, sigma, N)`, the same single
injection point as one-asset. Replace with our `income_process_Q.txt` loader
(≈ the 10 lines already in `swap_chain.py`); `nZ → 75`. No new discretization
work.

**What actually needs doing (on top of what we have):**

| task | effort |
|---|---|
| Swap our chain into `two_asset.make_grids`; retune `nB/nA` grids for `nZ = 75` | ~1 day |
| Re-calibrate: unknowns `(β, χ₁)`, targets `asset_mkt = 0` and liquid `B = B̄`; pick `B̄` and illiquid-wealth / `tot_wealth` targets (SCF liquid vs illiquid split — data we'd need to assemble) | ~2 days |
| The GE structure differs from ARS's: capital + investment (`Q, K, δ, εI`), sticky wages (`union`, `κw`), a bond/equity fund (`finance`, `ω`), distortionary labour tax. Understand it, decide which pieces to keep vs simplify toward KMV | ~2 days |
| Re-derive the direct/indirect **and** channel decomposition (§S3) for this model — it has *more* channels (capital return, investment), which is closer to KMV's Table 2 | ~2–3 days |
| Re-validate against KMV: steady-state wealth distribution (liquid + illiquid), MPC distribution, monetary decomposition | ~2 days |
| Re-run for USA/FRA/GER, check all three calibrate | ~1 day |

**Estimate: ~2 weeks of focused work**, most of it calibration and
decomposition, not model-building. Risk is low — the model and solver exist and
run; the uncertainty is in the wealth-split calibration data and in how much of
the richer GE structure to carry.

**Why not KMV's own replication code:** it is continuous-time
(HJB/KFE finite-difference, MATLAB), a different paradigm. Our discrete Markov
chain does not drop in — KMV parameterise the jump-drift process directly. The
one appealing version is feeding KMV our estimated `(λ, β, σ)` *directly*
(zero discretization error, KMV's validated solver), but that means a MATLAB
toolchain, their two-asset-only model, and re-doing the cross-country plumbing
in an unfamiliar codebase. Higher risk, ~3 weeks, only worth it if the
discretization (§8d) proves to be a first-order problem.

**Recommendation to put to the supervisor:** if the fork goes to two-asset, use
`sequence_jacobian.examples.two_asset` as the base (~2 weeks), not KMV's MATLAB.

---

## Open (step 2)

- **The one-asset calibration with our chains (§S4–S5, resolved for FRA/GER).**
  Adding KMV stochastic death to ARS's β-heterogeneity **makes FRA and GER
  calibrate** in one-asset (clean, MPC 0.20, IRF within 7% of baseline, KMV
  decomposition reproduced). **USA does not** even with β₂ pinned + death — its
  GRID chain lacks KMV's leptokurtic tail (`kurt Δ1y` 12.8 vs 17.8). Proposed
  resolution: FRA/GER on one-asset + death; USA validated by the S1 baseline
  (KMV's own chain), with the GRID-USA chain reported as a sensitivity. Poland
  is the residual risk (needs two-asset only if the Polish process is
  USA-shaped). Supervisor question, but the scope is now small.
- ~~Decomposition identity has a ~10% residual with the death block~~ **FIXED**
  (§S5, `death.py` shockable-transition MRO) — identity now holds to 1e-6.
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
