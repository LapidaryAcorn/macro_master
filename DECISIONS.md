# Modelling decisions and their justifications

Running log of choices made in the project, with the reasoning behind each, so
they can be written up in the thesis rather than reconstructed from memory.

Format per entry: **what was decided → why → how it was verified → what it
implies for the write-up.**

---

## 0. Goal (revised 2026-10-11): process *shape*, not Poland's decomposition alone

**Framing.** Supersedes §0a below. The thesis no longer asks "what does the
Poland decomposition look like" with USA/FRA/GER as method validation. It asks:
**how does the shape of the income process affect HANK outcomes, holding
everything else fixed?** At a common, shared set of calibration targets (the
same `A = 20`, SCF Lorenz curve, aggregate-MPC target used for USA/FRA/GER in
step 2), swap in the Polish estimated process and a reference process (US, or
whichever validates cleanly) into the *same* model and calibration procedure,
and compare MPC, the hand-to-mouth share, and monetary/fiscal transmission
across the two. **The object of interest is the difference the process shape
induces, not Poland's level in isolation.**

**Why this sidesteps the missing-Polish-target problem.** The open item at the
bottom of the old decisions list ("Poland's wealth moments are not in hand")
was a blocker *only* if the plan was to calibrate Poland's own `beta` to
Poland-specific wealth data. Under the comparison framing, every process is
run through the *same* calibration targets — only the income process differs
between runs — so no Poland-specific calibration target is needed for the
headline result. (A Polish-specific calibration remains a natural *robustness*
check, not a precondition.)

**Relation to work already done.** STEP2_LOG §S4–S5 (the chain swaps: USA/FRA/
GER processes dropped into the ARS one-asset model, re-calibrated to the same
targets, MPC/wealth/IRFs compared) is *already* a pilot run of this exact
exercise — it was framed there as "validation", but it is structurally
identical to the comparison the thesis now makes central. Worth stating
explicitly in the write-up: the comparison was rehearsed once, before Poland
existed, and Poland is now a third point on the same exercise.

**Sections whose stated rationale leaned on the superseded framing (§0a) and
should be reworded, or at least reread with this in mind, before the final
write-up:**
- **§0a** itself — superseded; kept below for the record, not deleted, because
  the sub-flags it raises on §8/§11/§12 are still partly live.
- **§8** ("re-run all three countries for consistency of method") — the same
  action is now justified *more* strongly: an identical pipeline across
  countries is a precondition for attributing HANK differences to process
  *shape* rather than to pipeline inconsistency.
- **§11** (common 2001–2016 window) — under a shape comparison a common window
  arguably matters *more*, not less (a confound: differences in the macro
  period each country's process was estimated over could themselves drive HANK
  differences, confounding "shape"). **Tension, recorded, not resolved:** Task 1
  (2026-10-11) uses Poland's full 2000–2020 GRID window rather than restricting
  to 2001–2016, because (i) the point of Task 1 is Poland's own best estimate,
  not a window-matched one, and (ii) the common-window machinery was built for
  a different purpose (cross-country *validation* comparability). If the
  supervisor wants the shape comparison itself to be window-matched, re-estimate
  Poland (and possibly USA/FRA/GER) on the shared window as a robustness check —
  flagged in the open-decisions list, not done by default.
- **§12** (GRID-USA ≠ KMV's published column) — reread "validity standard is
  consistency of method" as "the US estimate now *feeds* the comparison,
  it does not merely validate the pipeline in isolation"; the substantive
  finding (GRID-USA's shape fails to support the KMV mechanism at all, §S4–S5)
  becomes a *result*, not a caveat on a validation exercise.
- Anywhere the phrase "Poland is the thesis target, USA/FRA/GER validate the
  method" appears (in this file or STEP2_LOG.md), read it as "each country's
  estimated process is one point in a comparison whose point *is* the
  comparison" — including Poland.

**New open item this creates.** For the comparison to isolate *shape* rather
than *scale*, the countries being compared should probably sit at a common
ergodic dispersion (as §8b's `ErgodicVarConstraint` already does for GER, and
as USA now also uses) — otherwise "Poland has a different MPC" could be a
dispersion-scale artifact, not a shape one. Whether to apply the same
ergodic-variance pin to Poland by default, or to run both a free and a pinned
Polish estimate and compare, is a supervisor question; raised again under Task
1 below once the Polish estimate is in hand.

**For the write-up.** Open the thesis with this framing, not the
Poland-decomposition one. State plainly that the design routes around the
missing-Polish-wealth-target problem by comparing *within* a shared
calibration, rather than calibrating Poland in isolation against data that do
not exist yet.

---

## 0a. [SUPERSEDED 2026-10-11 — see §0] Scope: the thesis target is Poland; USA/FRA/GER are method validation

**Framing (set 2026-08-30).** The thesis studies **Poland only**. The USA, France
and Germany results are not findings in their own right — they exist to
**validate the pipeline**: evidence that it reproduces KMV on US data and behaves
sensibly on other countries before it is applied to Polish data (which are not
yet in hand). There is no cross-country comparison as a research question.

**Sections written under the old (cross-country) framing** — their *decisions*
stand, but their stated *rationale* leans on cross-country comparability and
should be reworded once the Polish data arrive and the validation framing is
final:

- **§8** ("re-run all three countries", not just the two with a binding bound) —
  justified by "consistency of method across countries". Under the new framing
  the same action is justified more simply: all three validation runs should go
  through the identical final pipeline.
- **§11** (common 2001–2016 window) — entire rationale is cross-country
  comparability. As a validation exercise the window choice matters much less;
  keep it for tidiness or note it is immaterial.
- **§12** (GRID-USA ≠ KMV's published US column; "validity standard is
  consistency of method across countries") — the validity standard is now
  "reproduces KMV's *mechanism* and moments on comparable US data", not
  cross-country consistency. The sample-difference explanation is unchanged.

Do not rewrite these yet — flagged here so the rewording is not forgotten.

---

## 1. Simulation design: additive jumps, not resets

**Decision.** Jumps add to the existing component (`z_j → z_j·e^{-β} + ε`),
rather than replacing it (`z_j → ε`).

**Why.** This is KMV's specification (eq. 30–31): `dz_j = -β_j z_j dt + ε_j dN_j`.
The `dN_j` term is additive by construction.

**Verification.** Tested both variants at KMV's published Table 3 parameters.
Only the additive version reproduces their reported model moments; the reset
version misses badly (e.g. var(Δ1y) ≈ 0.40 against a target of 0.23).

**For the write-up.** Not a modelling choice so much as a correctness
requirement — worth one sentence noting the distinction, since a reset
specification is a natural misreading of "jump".

---

## 2. Lifecycle panel, not the ergodic distribution

**Decision.** Workers enter at `z = 0` and are followed for 36 years
(≈ ages 25–60). Moments are computed on this panel, not on draws from the
process's stationary distribution.

**Why.** With `λ₂ = 0.007` (quarterly) the persistent shock arrives roughly
once every 38 years — about once per working life. A worker observed over a
career has therefore experienced only a handful of persistent shocks, whereas
the ergodic distribution represents a population that has experienced
arbitrarily many. The two are not interchangeable when the shock is that rare.

**Verification.** The ergodic version gives var(log y) ≈ 0.9 against KMV's 0.70;
the lifecycle panel gives ≈ 0.69. This also matches the empirical sample the
moments come from (GKOS follow prime-age workers over careers, not a synthetic
stationary population).

**For the write-up.** This is a substantive point about matching the simulation
design to the sampling design of the data. It also explains a set of numbers
that look odd in the output: the exported chains' *ergodic* var(log e) is far
above the 36-year *panel* var(log y) they were fitted to — USA ≈ 1.79 vs 0.96,
FRA ≈ 1.04 vs 0.48, GER ≈ 1.04 vs 0.55. Same process, different stage of
convergence, not an inconsistency. It also means step 2, if solved on the
stationary distribution, confronts households with markedly more earnings
dispersion than the estimation targets show — see §8b and the USA note in the
open-decisions list.

---

## 3. High-frequency simulation with time aggregation

**Decision.** Simulate at 6 sub-steps per quarter, accumulate the earnings
*flow* `exp(z)` within each year, and take logs of the annual total.

**Why.** The data moments are on *annual* earnings, which are the time
aggregate of a flow, not a snapshot of `z` at one instant. Simulating at a
coarse step and reading off a point-in-time value misrepresents both the
variance and the tail behaviour of annual earnings.

**Verification.** Convergence check across step sizes: 1 step/quarter gives
kurt(Δ1y) ≈ 12.7, 3 steps ≈ 15.0, 6 steps ≈ 14.6, 13 steps ≈ 14.7 — converged
by 6, which is the default.

**For the write-up.** Standard time-aggregation point, but worth stating
explicitly because it materially changes the fitted moments.

---

## 4. Objective function: squared *relative* deviations

**Decision.** Minimise `Σ ((simulated − target)/target)²` over the 8 moments,
with equal weights.

**Why.** The moments live on very different scales: variances ≈ 0.2, kurtoses
≈ 15. Under absolute deviations the kurtosis terms would dominate the objective
purely by magnitude, and the optimiser would effectively ignore the variances.
Dividing by the target puts every moment on a "percent off" footing.

**Alternative not taken.** A proper SMM weighting matrix (e.g. the inverse of
the moment covariance matrix) would be the textbook choice and would also
downweight imprecisely-estimated moments. It was not used because GRID does not
publish standard errors for the reported statistics, so the weighting matrix
cannot be constructed from the available data.

**For the write-up.** State the objective and note the weighting-matrix caveat
— an examiner familiar with SMM will look for it.

---

## 5. Log-scale parameter search

**Decision.** The optimiser searches over `log θ`, not `θ`.

**Why.** The parameters span two orders of magnitude (`λ₂ ≈ 0.007` vs
`β₁ ≈ 0.76`). In levels, a search step large enough to matter for `β₁` would
step clean over the entire plausible range of `λ₂`. In logs, a step means "a
given proportional change", which is scale-appropriate for every parameter.

**For the write-up.** One line in the estimation section.

---

## 6. Common random numbers in the simulator

**Decision.** Fixed seed and fixed RNG-call structure, so every objective
evaluation uses the same underlying random draws.

**Why.** Without this, two evaluations of the *same* parameters return
different objective values, and the optimiser cannot tell a genuine improvement
from a lucky draw — it chases simulation noise. With CRN the objective is a
deterministic, reasonably smooth function of the parameters.

**For the write-up.** Standard SMM practice; one sentence, but it is the reason
the optimisation converges cleanly at modest simulation sizes.

---

## 7. Parameter bounds separate the two components

**Decision.** `λ₁ ∈ (0.01, 0.6)`, `β₁ ∈ (0.05, 8.0)` (originally 4.0, see §8),
`λ₂ ∈ (0.001, 0.10)`, `β₂ ∈ (0.0001, 0.20)` (lower bound reduced from 0.001,
see Note), both `σ ∈ (0.20, 3.5)`.

**Why.** Three roles. (i) Technical: differential evolution requires a bounded
box. (ii) Identification: the near-non-overlapping ranges for `λ` and `β`
enforce the intended interpretation — component 1 frequent and fast-decaying,
component 2 rare and near-permanent. Without them the optimiser can relabel the
components or collapse both toward the middle, fitting the moments while
destroying the economic content. (iii) Validity: rates must be strictly
positive (`β = 0` makes the stationary variance `λσ²/(2β)` diverge).

**Note.** For GER, `β₂` is calibrated rather than estimated (see §8b), so its
bound does not apply there. The `β₂` lower bound was reduced from 0.001 to
0.0001 so a weakly-identified near-permanent component is not clipped; in the
event the estimated `β₂` for USA and FRA (≈ 0.005–0.007) sits well inside the
original range and GER's is calibrated, so the change is precautionary. The two
`σ`s deliberately share a range. The components are distinguished by their
*rates*, not their shock sizes — KMV's own estimates have similar `σ`s (1.74 and
1.53), so identifying restrictions are carried by `λ` and `β`.

**For the write-up.** Report the bounds and state the check applied to them:
no estimate should sit *at* a bound (see §8).

---

## 8. Widening the `β₁` bound from 4.0 to 8.0, and re-running all three countries

**Decision.** After the full run, `β₁` hit its upper bound of 4.0 for FRA and
GER. The bound was widened to 8.0 and **all three countries** re-estimated, not
just the two affected.

**Why widen rather than accept.** A parameter resting exactly on a bound means
the optimiser wanted to go further and was prevented — the reported value is an
artefact of the constraint, not an estimate. An objective sweep showed a genuine
interior optimum near `β₁ ≈ 4–4.5`, confirming the bound was binding
artificially rather than protecting against a runaway.

**Why it happened for FRA and GER but not USA.** These countries have much lower
variance of annual changes (FRA 0.207, GER 0.148, against USA 0.326). To
reproduce year-to-year earnings that stable while still admitting transitory
shocks, the transitory component must decay very fast — i.e. high `β₁`. The USA,
with more volatile earnings, needed less: `β₁ ≈ 2.76`, comfortably interior.
So the binding bound is a symptom of a real cross-country difference in
earnings stability, not a numerical accident.

**Why re-run USA too, despite it being clean.** The thesis's validity standard
is *consistency of method across countries* (see §11). Passing one country
through a different pipeline than the other two — different bounds, and a
separately-fixed polish step — would undercut exactly the argument the
cross-country comparison rests on. The extra ≈ 50 minutes of computation buys a
clean methodological story.

**Related fix.** The same run exposed a bug in which the Nelder-Mead polish did
not respect the bounds, which is how a parameter could reach the boundary in the
first place. Fixed alongside.

**Decision taken:** widen to 8.0 and re-run all three countries.

**Resolved after the re-run (serious budgets, 2026-08-30).** `β₁` settled at an
interior optimum for both affected countries — FRA `β₁ = 4.27` (log-position
0.88 in the widened box), GER `β₁ = 5.34` (0.92) — and did *not* run to the new
bound of 8.0. An objective sweep holding the other parameters at their estimates
confirms a genuine interior minimum in each case: for GER, `f = 0.043` at
`β₁ = 4.5`, `0.035` at `5.3`, `0.042` at `6.5`, `0.057` at `8.0`. So `β₁` *is*
separately identified from the annual moments and no parameter needs to be
fixed; the earlier bound was simply too tight. USA is unchanged (`β₁ = 2.75`).
Differential evolution returned essentially the same estimates in the two
independent serious runs (USA to four decimals; FRA `β₁` moved 3.96 → 4.27 only
because the bound was released), which is the stability one wants to see. The
"fix `β₁`" contingency is therefore not triggered.

**For the write-up.** This is a good methodological passage: it shows the bound
check was applied, that a binding bound was diagnosed rather than ignored, and
that the cross-country pattern has an economic reading.

---

## 8b. Germany: `β₂` is not identified — calibrated, not estimated

> **Revision note.** This section was first written as "fix `β₂ = 0` (random
> walk)" on the basis of an earlier, coarser run in which both optimisers drove
> `β₂` to zero. A finer objective sweep then showed a *weak interior* optimum at
> `β₂ ≈ 0.0007` rather than a corner solution at zero, and — more importantly —
> exposed a consequence that a corner solution had hidden. The decision below
> supersedes that earlier one. The original reasoning about the German result
> being economically meaningful still stands; what changed is how the parameter
> is pinned down.

**Decision.** For GER, `β₂` is **calibrated, not estimated**. It is not fixed at
a number but tied to the other persistent-component parameters by a constraint:
`β₂ = λ₂ σ₂² / (2 V₂*)` with `V₂* = 0.9103`, KMV's ergodic variance of the
persistent component (`0.007 · 1.53² / (2 · 0.009)`). So `β₂` is recomputed from
the current `λ₂`, `σ₂` at every objective evaluation and the ergodic sd of the
persistent component is held at exactly `√V₂* = 0.954` regardless of where the
search moves. The remaining five parameters are estimated; USA and FRA keep the
full six-parameter specification. Implemented as `ErgodicVarConstraint` in
`estimate.py`, passed as `derive={"beta2": ...}`. At the GER estimate the
constraint yields **`β₂ ≈ 0.0066`** (`λ₂ = 0.0102`, `σ₂ = 1.08`) — an earlier
note said `≈ 0.005`, from a coarse sweep before the constrained re-run.

**Why `β₂` cannot simply be estimated here.** The German panel does not identify
it. The objective has only a weak interior optimum at `β₂ ≈ 0.0007`, i.e. it is
close to flat in that direction: a wide range of values fits the eight moments
about equally well. This is the §9 identification problem in its sharpest form —
with a half-life of ~250 years at the optimum, essentially no persistent shock
ever visibly decays within a 36-year panel, so the data carry almost no
information about the decay rate.

**Why the unconstrained estimate is unusable despite fitting well.** The
stationary variance of a jump-drift component is `λσ²/(2β)`, so as `β₂ → 0` the
implied ergodic dispersion explodes. At `β₂ ≈ 0.0007` the persistent component
has an ergodic **sd ≈ 2.5** in logs — an implied long-run earnings distribution
spanning several orders of magnitude, which is not a credible description of any
labour market. This is invisible in step 1, because the moments are computed on
a 36-year lifecycle panel that never approaches the ergodic distribution (§2).
It would not be invisible in step 2, which *is* solved on the stationary
distribution (confirmed — see weakness 2 below): that dispersion would drive the
wealth distribution, the MPC distribution, and possibly the solver's
convergence. The same mechanism broke the analytic variance correction (§8c).

**Cost of the restriction (measured on the re-run).** OOS objective goes from
`f ≈ 0.003` (USA/FRA territory) to `f ≈ 0.035` for GER — better than the
`≈ 0.08` first feared, but still an order of magnitude worse than the two
unrestricted countries. The deterioration lands almost entirely on **var(log
annual earnings)**: the discretized/continuous chain reproduces `0.55` against a
data target of `0.645`, a **15% underfit**. The 5-year change variance is `+10%`
(0.315 vs 0.285). Crucially the **kurtosis is unharmed** — `kurt Δ1y = 17.9` vs
target `17.8`, `kurt Δ5y = 11.2` vs `11.2` — so the leptokurtic feature the KMV
specification exists to capture (§14) survives the restriction intact. The
reading: pinning `β₂` costs the model some of Germany's *level* of earnings
dispersion but none of its *shape*.

**Honest framing: this is calibration anchored on the US.** Fixing the ergodic
sd to KMV's ≈ 0.95 anchors a German parameter to a US estimate — one level of
abstraction up from simply setting `β₂ = 0.009` (KMV's own value), but the
anchor is the same. The defensible argument is *not* "the data say so"; it is:
the German data do not identify `β₂`, so some value must be imposed from
outside; given that, it is better to impose it on an economically interpretable
quantity (long-run earnings dispersion) than on a parameter whose units carry no
intuition. In the thesis this belongs in the **calibration** section, not the
estimation section, and the cross-country anchoring should be stated plainly
rather than buried.

**Two weaknesses of this choice, recorded deliberately.**
1. *The anchor value is arguably wrong for Germany.* An ergodic sd of 0.95 comes
   from US data, and Germany has visibly lower dispersion on every observed
   moment (var log y 0.645 vs 0.956; var Δ1y 0.148 vs 0.326). A criterion
   scaled to Germany's *own* observed dispersion would be equally simple and
   less arbitrary. Worth revisiting.
2. *Could the problem be an artefact of the step-2 solution concept?* **Checked
   (2026-08-30) — no, the restriction is genuinely needed.** The Auclert–Rognlie–
   Straub `annual-review` code solves an **infinitely-lived household block on
   its ergodic distribution** (`household.py`: `sj.utilities.discretize.stationary(Pi_e)`,
   no lifecycle / OLG structure). So the ergodic distribution *is* what
   households face, and the unconstrained `β₂ ≈ 0.0007` (ergodic sd ≈ 2.5) would
   materialise in the model. The faithful estimate cannot be kept; §8b stands.
   (Bonus finding: that code already loads a KMV chain from text files
   — `inputs/kmv_process/ymarkov_combined.txt` + `ygrid_combined.txt`, via
   `expm` then row-normalise then stationary dist then exp-grid normalised to
   mean 1 — so swapping in `income_process_*.txt` is close to drop-in. Note
   their file is a *generator*; our `_P.txt` is already `expm(Q)`.)

**What survives from the original reading — stated carefully.** The unconstrained
German fit prefers a **weaker** decay of the persistent component than the USA
or France: its objective keeps improving as `β₂` falls to ~0.0007 (half-life
~250 years), whereas USA and France have genuine interior optima at
`β₂ ≈ 0.005–0.007` (half-life ~25–35 years). So there *is* a cross-country
ordering — German persistent shocks look closer to a pure random walk. But do
not overstate it: 25–35 years is already "permanent over a working life", so the
economically relevant statement is "all three persistent components are
near-permanent, Germany's most of all", not "Germany is a random walk and the
others mean-revert". And under the calibrated spec the ordering is removed by
construction (all three pinned to the same ergodic dispersion). A
near-random-walk permanent component is also the canonical specification in the
permanent–transitory literature (Blundell, Pistaferri & Preston; MaCurdy), where
it is typically assumed rather than estimated. An institutional reading —
sectoral bargaining, employment protection, lower job mobility — is
interpretation to discuss with the supervisor, not something the estimation
establishes.

**For the write-up.** Frame as: *the same estimation procedure is applied to all
countries; for Germany it returns a persistent component whose decay rate the
panel cannot identify, and which drives the implied long-run distribution to
implausible values, so `β₂` is calibrated rather than estimated.* Report the fit
cost and the anchoring choice explicitly.

**Direct confirmation added 2026-10-12:** §21 shows `β₂` cannot be recovered
even from noise-free synthetic data generated from a known truth, with the
full estimation budget and a neutral start. That is the strongest available
evidence for calibrating rather than estimating it here.

---

## 8c. Chain export: simulation-based scale calibration (`match_var_log=True`)

**Decision.** Enable `match_var_log=True` when exporting the discretized chain,
for all three countries.

**Why it became necessary.** The analytic variance correction of §13 rescales
each component by `√(v_theory/v_chain)` with `v_theory = λσ²/(2β)`. As `β₂ → 0`
(§8b) that denominator vanishes and `v_theory` explodes — in the limit a random
walk has **no stationary distribution** at all, so the quantity the correction
targets does not exist. The symptom, seen on the first (free-`β₂`) run for GER,
was a chain with var(log y) ≈ 1.41 against a continuous benchmark of ≈ 0.64.
Constraining `β₂` (§8b) pulls the ergodic variance back into a finite range and
removes the blow-up, but the general point below still stands: any correction
anchored on the ergodic variance is fragile as `β₂` gets small.

**Why this is a domain-of-validity issue, not an implementation bug.** No
refinement of the analytic scheme can fix it: any correction anchored on the
*ergodic* variance is the wrong instrument for a process that has no ergodic
distribution. The simulation-based alternative calibrates one global scale
factor against the variance of log annual earnings on the 36-year lifecycle
panel — a quantity that is well-defined for both mean-reverting and random-walk
processes.

**Why enable it for USA and FRA too, where the analytic method works.** Two
reasons. First, it is the better method there as well: §15 already documents
that the analytic correction leaves discretized variances ≈ 7% below the
continuous process, precisely because it targets the ergodic variance while the
moments come from a finite lifecycle panel — which is exactly what
`match_var_log` measures instead. Second, this step produces the files that feed
step 2, so keeping the export procedure identical across countries matters for
the same reason as §8 and §11.

**Cost.** Two extra simulations per country. Negligible.

**Verified after the run (2026-08-30).**
(i) *Discretized var(log y) now tracks the continuous process for every
country:* USA `0.955` (continuous `0.953`), FRA `0.479` (`0.476`), GER `0.551`
(`0.551`) — all within ~1%, against gaps of several percent under the analytic
correction alone. Because `match_var_log` targets the continuous model's
lifecycle-panel variance, and for these estimates that variance ≈ the data
target, the discretized column now also lands on the *data* — except GER, which
lands on `0.55` (15% below the `0.645` data target) for the separate reason in
§8b, not because of this step.
(ii) *Kurtosis is untouched by the scale factor,* as predicted. The run with
`match_var_log=True` gives USA `kurt Δ1y = 14.36` in the discretized column
against `14.33` in the earlier run without it — identical to rounding. The
global scale factor is a pure rescaling and kurtosis is scale-invariant;
confirmed empirically, not just assumed. (Grid *extent* does move kurtosis —
75 → 351 states pushes USA `kurt Δ1y` from 14.4 to 17.4 — but that is §14
territory, unrelated to the scale factor.)
(iii) *GER's chain variance came down from `1.41` to `0.55`.* Two things
contributed: constraining `β₂` (§8b) removed the exploding `v_theory`, and
`match_var_log` then pinned the panel variance to the continuous `0.55`. The
`match_var_log` step remains the general method — it would still be needed if
`β₂` were ever left free.

**For the write-up.** A clean methodological point: the analytic correction has
a domain of validity (mean-reverting processes with finite stationary variance),
the German estimate falls outside it, and the simulation-based calibration is
the general method covering both cases.

---

## 8d. The discretized chains understate the 1-year change variance (~30%)

**Observation (surfaced by the re-run).** With `match_var_log=True` the exported
chains reproduce var(log annual earnings) and the 5-year change dynamics well,
but the **1-year** change variance comes out about 30% low, with a compensating
rise in its kurtosis:

    var Δ1y     data    continuous   discretized chain
    USA         0.326     0.327          0.225
    FRA         0.207     0.206          0.122
    GER         0.148     0.142          0.085

    kurt Δ1y    data    continuous   discretized chain
    USA         12.8      12.1           14.4
    FRA         15.2      15.3           21.2
    GER         17.8      17.9           25.2

**Not a grid-resolution problem — a non-obvious result worth stating.** Refining
the grid does not help and slightly hurts: going 75 → 189 → 351 states leaves
USA var Δ1y stuck at ≈ 0.23 while kurt Δ1y drifts from 14.4 up past 17,
overshooting the continuous 12.1. More points widen the represented tails
without restoring the missing central variance. So the usual lever for
discretization error (a finer grid, §14) does nothing here; the problem is on
the time axis, not the state axis.

**Diagnosis — a frequency mismatch.** The estimated transitory decay rates are
high: `β₁ ≈ 2.75 / 4.27 / 5.34` per quarter, a half-life of a few weeks. The
continuous simulation runs at 6 sub-steps per quarter and time-aggregates the
flow, so it captures the within-quarter path of the transitory component that
contributes to annual-change variance. The exported chain is a *quarterly*
transition matrix (`expm(Q)` over one quarter); a shock that mostly decays
inside a quarter is barely resolved at that sampling frequency, so its
contribution to the annual change is lost. This is intrinsic to handing a
quarterly Markov chain to a quarterly model when part of the income process
lives at a higher frequency.

**KMV do not have this problem — checked, and it confirms the diagnosis.**
Running the same discretization at KMV's own Table 3 parameters (`β₁ = 0.761`, a
transitory half-life near one quarter):

    var Δ1y at KMV params    continuous   discretized (analytic)   discretized (match_var_log)
                               0.228          0.199  (−13%)              0.216  (−5%)

Both bracket KMV's published Table 2 discretized column (0.21 against data 0.23,
−9%). So the method itself is faithful — it reproduces KMV's own discretized
moment. The ~30% gap for our countries is not a flaw in the method; it is a
consequence of the parameter regime our data select. **This deserves an explicit
sentence in the thesis:** the GRID data for all three countries force a
transitory component that decays far faster than KMV's (a symptom, like the
`β₁` bound in §8, of much lower annual-change variance — FRA 0.21, GER 0.15 vs
KMV/USA ~0.33), and a quarterly chain cannot carry a component that fast. The
discretization consequence follows from the empirical finding, not the other way
round.

**One phenomenon, two consequences — the §8 connection.** High `β₁` is a single
empirical fact about the GRID data, and it shows up twice in this log: in §8 it
pushed the estimate against the parameter bound (bound raised 4.0 → 8.0, genuine
interior optimum at 4.3–5.3), and here in §8d it breaks the quarterly
discretization of the transitory component. Both trace to the same cause —
year-to-year earnings in France and Germany are much more stable than in KMV's
US sample, which the model can only reproduce with a very fast-decaying
transitory shock. Worth writing up as one point with two symptoms rather than
two separate curiosities.

**Why it matters for step 2.** KMV's monetary-transmission results run through
precautionary saving and the MPC distribution, which respond to earnings risk at
roughly the annual (business-cycle) horizon. A chain carrying only ~70% of the
1-year change variance will understate that risk and bias the step-2 results in
a direction that is hard to sign a priori.

**Decision (to implement in step 2).** Discretize component 1 **directly as a
near-iid quarterly shock** — its own small Markov chain built from the quarterly
distribution of `z₁`, rather than routed through the drift-diffusion generator.
At `β₁ ≈ 3–5` the transitory component is almost iid across quarters anyway
(>90% decay per quarter), so a direct construction should carry its variance
without the frequency loss. The Kronecker structure of the combined chain is
unchanged. **Verify against the continuous var Δ1y after building it**, and
confirm var(log y) and the kurtoses do not regress.

*Fallback if that does not close the gap:* add a separate iid transitory shock
inside the HANK model on top of the persistent chain, calibrated to the var Δ1y
target — close to KMV's own structure, but it puts one income parameter in the
model layer rather than the estimated process, so it is second choice.

*Not chosen:* simply accepting the ~30% shortfall. The missing variance is at
the 1-year horizon, which is exactly where precautionary saving responds; it is
not innocuous high-frequency noise.

**For the write-up.** Belongs in the discretization-limitations subsection (§15),
with the resolution recorded in the step-2 methodology.

---

## 9. Two optimisers as an identification check, not a contest

**Decision.** Every country is estimated with both differential evolution and
TikTak, and the two results are compared — specifically on the persistent
component (`λ₂`, `β₂`, `σ₂`).

**Bottom line (after the re-run — read this before writing §9 up).** The
"two independent optimisers agree" form of the check **does not hold** for this
objective. TikTak's sequential Nelder-Mead local phase systematically stalls
short of the optimum (piecewise-flat CRN objective, below), landing at
`β₁ ≈ 2.2` for FRA/GER while DE reaches `β₁ ≈ 4.3 / 5.3` with a strictly better
fit — and it still stalls with a warm start and `n_sobol = 2048`. So TikTak is
**not** a reliable independent check here and must not be presented as one.

The identification evidence that *does* stand is:
  (a) a **line scan** between the DE and TikTak solutions showing the objective
      declines monotonically from one to the other — a single valley, not two
      basins, so this is a convergence gap, not a multi-valley identification
      problem; and
  (b) **DE's stability across two independent serious runs** (same estimate to
      3–4 digits), plus interior objective sweeps in every parameter (§8).

Write it up that way. Do not claim verification by optimiser agreement.

**Why.** The persistent component is weakly identified: with `λ₂ ≈ 0.007` a
36-year panel contains few persistent shocks, so several parameter
configurations fit the 8 moments almost equally well. Two structurally
different search algorithms landing in the same place is evidence of a genuinely
well-identified optimum; landing in different places is a flag to investigate —
though, as the re-run showed, "different places" can still mean "same valley,
one optimiser stalled" rather than a genuine multi-valley objective, so the
divergence has to be diagnosed (line scan between the two solutions), not just
noted.

**Evidence this matters.** On early short-budget runs, DE found a
near-random-walk persistent component (`λ₂ ≈ 0.047`, `β₂ ≈ 0`, `σ₂ ≈ 0.40`)
that fit the moments about as well as a KMV-like configuration but means
something economically quite different. See the re-run outcome below for how the
two optimisers compare at serious budgets.

**Important distinction.** Optimiser choice affects the *reliability of finding*
the best-fitting point. Identification affects whether that point is
*economically unique*. No optimiser can separate two observationally equivalent
parameter sets — that would need an additional identifying moment (e.g. a
longer-horizon autocovariance, which GRID reports).

**Compare the implied quantities, not just the raw parameters.** The parameters
are not all economically meaningful one-for-one — `λ₂`, `β₂`, `σ₂` matter mostly
through the combinations they imply (ergodic dispersion of the persistent
component `√(λ₂σ₂²/2β₂)`; the transitory half-life `ln2/β₁`). The agreement
check is on those. Reporting it purely as "`β₂` differs by 27%" both overstates
a benign disagreement and would let a genuine one hide behind offsetting
parameter moves.

**Outcome on the full re-run (serious budgets, 2026-08-30).**

*Persistent component — agrees, on the quantity that matters.* FRA and GER: the
two optimisers agree on the raw parameters (`λ₂`/`σ₂` within 9%). USA: `β₂`
differs by 27.5% (DE 0.0066, TikTak 0.0048) and `σ₂` by 13% (1.32 vs 1.14) —
but the **implied ergodic sd of the persistent component is 1.22 (DE) vs 1.35
(TikTak)**, i.e. the parameters trade off against each other and the
economically meaningful quantity moves only ~10%. On that reading USA passes the
check. DE additionally reproduced its own persistent-component estimate across
the two independent serious runs to three digits.

*Transitory component — the two optimisers land apart, but on the same valley.*
`β₁` is 4.27 (DE) vs 2.01 (TikTak) for FRA and 5.34 vs 2.42 for GER, with `λ₁`
differing too, and DE fits 3–4× better out of sample. A better OOS fit alone is
*not* enough to declare that resolved — that is the reasoning the §9 check
exists to resist — so it was tested two ways.

1. *Re-run TikTak with symmetric information.* `scripts/tiktak_recheck.py` runs
   TikTak again for FRA and GER with a KMV warm start (same as DE) and
   `n_sobol = 2048`, `n_local = 15`. It **still does not reach the high-`β₁`
   solution** — it lands at `β₁ ≈ 2.2 / 2.3`, OOS `f = 0.0045 / 0.0397` (barely
   better than before, still well above DE's `0.0013 / 0.0353`). So the gap is
   not fixed by budget or by the warm start.

2. *Scan the objective on the straight line between the two solutions.* In
   log-parameter space, from TikTak's point to DE's point, the objective is
   **monotonically decreasing** for both countries, with no intermediate hump:

       FRA:  f = 0.0042 → 0.0036 → 0.0032 → 0.0027 → 0.0022 → 0.0017 → 0.0013
       GER:  f = 0.0393 → 0.0388 → 0.0387 → 0.0383 → 0.0378 → 0.0369 → 0.0350

   (FRA strictly monotone; GER monotone bar one 0.0001 stair-step near the flat
   top — i.e. within simulation noise, no real hump.) There is a clean downhill
   path connecting the two solutions. They are **not separate basins** — DE's is
   simply the better-converged point on a single valley, and TikTak stalled
   partway down it. Data in `output/tiktak_line_scan.json`.

*Why TikTak stalls.* The SMM objective is piecewise-flat: with a fixed 50k-worker
CRN panel, the simulated moments change in discrete jumps as parameters cross
thresholds that move a worker across a grid bin or add/drop a jump. The DE
iteration log shows this directly — `f` sits constant for 5–10 generations at a
time. Sequential Nelder-Mead's simplex collapses on those flat patches and
reports convergence; DE's population steps across them and keeps descending.
This is a known Nelder-Mead failure mode on stair-shaped objectives, not a
statement about identification.

*Conclusion.* The transitory component **is** identified — a single valley
(line scan), a stable optimum across two independent DE runs, and interior
objective sweeps in every parameter (§8). DE is canonical. TikTak's
Nelder-Mead local search is not a reliable optimiser for this objective; the
meaningful cross-check here is the line scan plus the DE cross-run stability,
not TikTak's endpoint. For a genuinely independent second optimiser the local
phase would need a pattern-search method (or a much larger, smoother
simulation) — recorded in the open-decisions list.

**For the write-up.** State plainly: identification of the estimated process
rests on (a) the line scan showing a single valley connecting the two
optimisers' solutions and (b) DE's reproducibility across independent runs and
the parameter-by-parameter objective sweeps — **not** on the two optimisers
agreeing, because TikTak's local search is unreliable on this stair-shaped
objective. Report the persistent-component comparison in terms of the implied
ergodic dispersion, not the raw `β₂`. The separate question of whether the 8
moments *could* identify a genuinely different economic configuration (as
opposed to whether the optimiser finds the best one) is unchanged and would
still need a ninth moment — keep that distinction explicit.

**Amendment (2026-10-12, after §21).** The recovery exercise weakens one claim
made above. §9 concludes "DE is canonical" and rests the identification
argument partly on "DE's stability across two independent serious runs" — but
both of those runs were warm-started at the KMV Table-3 values on real data.
§21 shows that on GER-shaped data, DE started from that same neutral point
lands at `β₁ ≈ 1.84`, in the same low-`β₁` neighbourhood where TikTak's
Nelder-Mead stalled (`β₁ ≈ 2.2–2.4`), in 4 of 4 draws. So the honest statement
is narrower than §9's: DE is the better optimiser on this objective and
reproduces its own answer run-to-run, but it is not demonstrably able to find
the global optimum from an arbitrary start for a near-random-walk persistent
component. For GER specifically the canonical estimate is defended not by
DE's reliability but by a direct objective comparison against the alternate
basin on real data (§21: 2.6× better). What stands unchanged: the line scan
showing a single valley, and the parameter-by-parameter sweeps (§8).

**Line scan extended past the DE point (2026-10-12).** The open question this
amendment raised — is the canonical point actually the floor of the valley,
or does the objective keep falling past it in the same direction? —
`scripts/extend_line_scan.py` answers directly: continue the same
log-parameter-space line from the TikTak point through the DE point (`a=0` to
`a=1`, exactly reproducing `output/tiktak_line_scan.json`'s construction) out
to `a = 1.25, 1.5, 1.75, 2.0` — i.e. 25-100% further in the same direction,
on real targets, same CRN. (The a=0/a=1 reproduction came out close but not
bit-for-bit identical to the original file — FRA `f=0.00418` vs original
`0.00420` at `a=0`, `0.00121` vs `0.00127` at `a=1`; GER similarly ~1-2% off.

**Correction (2026-10-12): this is not environment drift — checked directly,
twice.** First guess was a numpy RNG change since the original scan predates
this session. Wrong. Reproduced FRA's `a=1` (DE) point independently on a
*different machine* (numpy 2.5.3) running the code as it stood at commit
`6c060bd`, and got `f=0.0012070` — matching this session's own
`extend_line_scan.py` value (`0.0012069597826309614`) to 7 significant
figures, and *not* matching the archived `0.001269` in
`tiktak_line_scan.json`. So today's code is perfectly reproducible
machine-to-machine and commit-to-commit; the archived value simply came from
whatever uncommitted, since-lost script originally produced
`tiktak_line_scan.json` (confirmed it was never checked in — there is no
script in the repo history that generates that file), running with settings
that are not quite the same as `extend_line_scan.py`'s. The ~1-2% gap is a
different-script artefact, not a reproducibility problem with current code.

**Separately, and worth stating on its own: Monte Carlo noise at this point
is substantial relative to `f` itself.** Evaluated the objective at FRA's DE
point under 5 different seeds and 3 different simulated-panel sizes (holding
parameters fixed): `f` ranged from `0.00121` to `0.00233` (mean `0.00160`,
sd `0.00036`, i.e. **sd/mean ≈ 23%** of the point's own objective value).
Practical implication: a ratio near 1 between two objective values at this
scale — e.g. §22's GER `V3` result (`0.0028` vs `0.0025`, a `1.1×` ratio) —
sits well inside this noise band and should not be read as a real
difference. A ratio of several-fold or more (§22's `V1`/`V2` results, both
`>100×`) is far outside it and is a real signal.)

**Result: the objective rises sharply and monotonically past the DE point for
both countries — the canonical estimate is the valley floor, not a point on a
still-descending slope.**

| | FRA `f` | GER `f` |
|---|---|---|
| a=1.00 (DE, canonical) | 0.00127 | 0.03500 |
| a=1.25 | 0.00264 (+108%) | 0.03893 (+11%) |
| a=1.50 | 0.00907 (+614%) | 0.05338 (+53%) |
| a=1.75 | 0.02486 (+1857%) | 0.08898 (+154%) |
| a=2.00 | 0.05531 (+4257%) | 0.15674 (+348%) |

FRA's objective is already more than double by `a=1.25` and rises over 40-fold
by `a=2.00`; GER rises more gradually at first but still clearly and
monotonically, quadrupling by `a=2.00`. **This settles the open question from
the amendment above: the canonical point is the valley floor along the
direction connecting the two optimisers — it is not merely "better than
TikTak", it is a genuine local minimum extending well past it.** Combined
with the parameter-by-parameter sweeps (§8) and this extended scan, three
independent pieces of evidence now agree the canonical estimate sits in a
true (if shallow, per §21's recovery evidence) basin, not on an open slope.
Data in `output/tiktak_line_scan_extended.json`.

---

## 10. GRID sample selection: male, 25–55, residualized changes

**Decision.** Male, prime age 25–55, residualized (age/time effects removed)
changes, statistics averaged over years.

**Why.** To match the sample KMV target (GKOS's SSA extract). Pooling genders
is not comparable: an "All genders" extract inflates var(log y) and the change
variances and depresses the small-change fractions, because it mixes in
participation margins and part-time patterns. Prime age avoids
school-to-work transitions at one end and retirement decisions at the other,
both of which distort the variance and especially the tails without saying
anything about the income *process*.

**Verification.** Switching from an "All genders" to a male-only extract was
tested directly; the all-genders version gave var(log y) ≈ 0.93.

**For the write-up.** State the selection explicitly; it is the first thing a
reader will want to know when comparing to KMV.

---

## 11. Common year window 2001–2016 across all countries

**Decision.** All countries' moments computed over 2001–2016, the period
covered by all three in the GRID extract (native coverage: FRA 1991–2016,
GER 2001–2016, USA 1998–2019).

**Why.** Comparing France measured over the 1990s–2010s against the USA measured
over the 2000s–2010s would let differences in *when* each country was observed
contaminate the cross-country comparison. A common window removes that channel
by construction. The argument is methodological rather than empirical: it makes
the comparison defensible regardless of how large the trend effects turn out to
be.

**Empirical note.** The window in fact changes very little (USA var_log
0.945 → 0.956, FRA 0.488 → 0.476, GER unchanged since its coverage already is
2001–2016). So the choice costs almost nothing in precision and is worth making
for the cleanliness of the argument.

**For the write-up.** State the window and the coverage table, and note that
the restriction is immaterial empirically — this pre-empts the objection that
the window was chosen to produce a result.

---

## 12. GRID-USA is not expected to reproduce KMV's published US column

**Decision.** Treat the discrepancy as a sample-definition difference to be
documented, not a bug to be fixed, and adopt *consistency of method across
countries* as the validity standard.

**Why.** KMV target GKOS's SSA extract over ≈ 1978–2013; GRID-USA here is
1998/2001–2019 with GRID's own harmonized construction. US earnings dispersion
rose over these decades, so a later window mechanically gives higher variances
and lower fractions of small changes.

**Evidence it is a sample difference, not an error.** All eight moments move in
the direction a "more volatile sample" implies, and the small-change fractions
diverge by amounts that shrink as the threshold widens (−0.15 at 10%, −0.12 at
20%, −0.05 at 50%) — exactly what a wider change distribution produces. A coding
error would not produce that coherent a pattern. A likely second contributor is
GKOS's minimum-earnings threshold, which trims volatile low earners and is not
replicated here.

**Step-2 finding (2026-08-31) — the GRID-USA chain does not support the KMV
mechanism, and that is the informative part.** In step 2 (STEP2_LOG §S4–S5) the
GRID-USA chain does not merely *produce different numbers* in the one-asset
HANK — it **fails to calibrate at all**: no discount-factor configuration (with
or without `β₂` pinned, with or without stochastic death, and after rescaling
the ergodic variance to KMV's) can jointly hit `A = 20`, MPC `= 0.20`, and the
SCF Lorenz curve. The aggregate MPC collapses to ≈ 0.03. **It is not the
variance** (ruled out by direct test — rescaling the ergodic variance to KMV's
1.07 or to FRA/GER's 1.04 does not help). GRID-USA's earnings *shape* departs
from KMV's on several axes that all reduce the constrained fraction: thin tail
(`kurt Δ1y` 12.8 vs KMV 17.8), a more persistent transitory component
(`β₁ = 2.93` vs FRA 4.27 / GER 5.34), and more frequent transitory shocks
(`λ₁ = 0.22`). Households face a lot of medium-persistence income risk they
build precautionary buffers against, so `< 1%` end up at the constraint.
**France and Germany — whose GRID `kurt Δ1y` (15.2, 17.8) and process shape are
KMV-like — do calibrate** (one-asset + death, MPC 0.20, KMV decomposition
reproduced). So this is a documented result: *the shape of the earnings
process, not just its variance, determines whether it supports the KMV
transmission mechanism, and the GRID-USA sample's shape does not.*

**For the write-up.** State the sample difference openly and early. Then, in
step 2, state the mechanism consequence explicitly: the US validation runs on
KMV's own chain (the S1 baseline reproduces their 20/80 split); the GRID-USA
chain is reported as a sensitivity whose *failure* to support the mechanism is
itself the finding.

---

## 13. Discretization: correcting the variance bias analytically

**Decision.** After building each component's transition matrix, rescale its
grid by `√(v_theory / v_chain)`, where `v_theory = λσ²/(2β)` is the exact
stationary variance.

**Why.** The upwind finite-difference treatment of the drift is numerically
diffusive: probability mass hops a whole grid point at a time instead of
drifting smoothly, which inflates dispersion. At KMV's parameters the raw chain
has ergodic variance ≈ 1.27 against a theoretical 1.07 — about 19% too high.

**Why a correction rather than a finer grid.** The bias is systematic, not a
resolution artefact: it is essentially unchanged from a 33-state to a 465-state
grid. Refining the grid does not remove it.

**Why this correction is safe.** Kurtosis is scale-invariant, so rescaling
fixes the variances without disturbing the tail behaviour. Verified: per-component
chain variances now match theory to four decimals.

**Domain of validity.** This correction assumes the component *has* a finite
stationary variance. It breaks down for a random-walk component (`β = 0`), where
`λσ²/(2β)` diverges — see §8c, where the simulation-based calibration replaces
it.

**For the write-up.** Worth a short methodological note — it is a real finding
about the discretization scheme, not just plumbing.

---

## 14. Discretization: default grid 5×15 = 75 states, not KMV's 3×11 = 33

**Decision.** Default to 75 states, with `n1`/`n2` exposed as arguments.

**Why.** Kurtosis requires grid points far enough out to represent rare large
jumps. Against a continuous benchmark of kurt(Δ1y) = 14.6: 33 states give 10.9,
75 give 13.9, 189 give 14.7, 275 give 14.8. Below ≈ 75 the tails are visibly
clipped; above ≈ 189 nothing is gained.

**Why it matters beyond the table.** The discretized chain — not the continuous
process — is what the HANK model consumes. Tail risk in earnings drives
precautionary saving and the MPC distribution, which are the mechanisms KMV's
monetary-transmission results run through. A chain that understates the tails
would quietly bias the step-2 results.

**Trade-off to revisit in step 2.** 75 states instead of the usual 7–11 will
slow the model solve. If it proves too slow, 33 states remain available and
still benefit from the variance correction (kurtosis 11.9 rather than 8.1).

**For the write-up.** Report the grid size and the convergence table; it
demonstrates the choice was made on evidence rather than convention.

---

## 15. Known residual limitations (documented, not hidden)

- **Discretized variances ≈ 7% below the continuous process** under the
  analytic correction alone. Addressed as of §8c by enabling
  `match_var_log=True` for all countries, which calibrates the scale against the
  lifecycle-panel variance directly.
- **`frac Δ1 < 10%` too high in the discretized column (≈ 0.65 vs 0.55).**
  Intrinsic to discretization: in many quarters the chain does not change state
  at all, piling mass at near-zero annual changes. No grid refinement removes
  it.
- **1-year change variance ≈ 30% low in the discretized column** (§8d), because
  the estimated transitory component decays faster than the quarterly chain can
  resolve. Not a grid-resolution issue; resolution to be chosen in step 2.
- **1-year change kurtosis comes in substantially *higher* in the
  discretized column than in the continuous model — the opposite direction
  from the variance gap above, and country-varying in size.** Measured
  directly from `table2.csv` (discretized vs continuous, not vs data):
  POL +18.0%, USA +22.0%, FRA +38.1%, GER +40.8%. Both this and the `var Δ1y`
  shortfall trace to the same root cause — a quarterly chain under-resolving
  a transitory component that mostly decays within a single quarter — but
  show up as opposite-signed distortions in variance versus kurtosis of the
  same moment. Not country-specific; present, with varying magnitude, in
  every country's discretized chain.
- **Residual kurtosis gap in validation** (14.5 simulated vs KMV's 16.5 at their
  own parameters), attributed to GKOS sample details not modelled here —
  principally the minimum-earnings threshold.
- **Fractions of small changes are interpolated**, not observed. GRID reports
  quantiles, not `P(|Δ| < c)`; each percentile column gives a point of the CDF
  and the fractions follow from `F(c) − F(−c)` by linear interpolation. This is
  interpolation (the thresholds lie inside the percentile range), never
  extrapolation, and is the intended use of GRID's quantile statistics.
  Accuracy depends on the fineness of the percentile grid, which is high here.
- **`β₂` is not recoverable for a near-random-walk persistent component**
  (§21): on synthetic data generated from a known truth, GER's `β₂` comes
  back ~82% low in every draw, and Poland's fails in 1 of 4. Reported as an
  identification limit of the 8-moment SMM design, not of the implementation.

**For the write-up.** A short "limitations" subsection. Stating these is
stronger than leaving them to be discovered.

---

## 16. Step-2 interface notes (from reading the Auclert `annual-review` code)

> Step 2 has its own running log: **`hank_step2/STEP2_LOG.md`**. So far:
> S1 baseline reproduced & verified; S2 frequency = quarterly (six ways); S3
> one-asset reproduces KMV's 20/80 split (composition differs — see below);
> **S4** the raw exported chains do not calibrate one-asset for USA/FRA
> (`dβ → 0`, MPC ≈ 0.03); GER works. **S5**: KMV's ergodic var is 1.07 not 0.85
> (FRA/GER at 1.04 already match); death *alone* does not generate MPC in
> one-asset — but **death + ARS β-heterogeneity together DO make FRA and GER
> calibrate** (clean, MPC 0.20, KMV decomposition reproduced). **USA** does not
> calibrate even with `β₂` pinned + death — its GRID chain lacks KMV's
> leptokurtic tail. Proposed: FRA/GER on one-asset + death; USA validated by the
> S1 baseline, GRID-USA a sensitivity. **S6** costs the two-asset fallback
> (~2 wk, SSJ `examples/two_asset`). The notes below are the pre-step-2 reading
> of the code and still stand.

Recorded 2026-08-30 from `household.py` in `shade-econ/annual-review`, so step 2
starts from facts rather than assumptions. The relevant block is:

    Pi_e = linalg.expm(np.loadtxt('inputs/kmv_process/ymarkov_combined.txt'))
    Pi_e /= np.sum(Pi_e, axis=1)[:, np.newaxis]
    pi_e = sj.utilities.discretize.stationary(Pi_e)
    e_grid_short = np.exp(np.loadtxt('inputs/kmv_process/ygrid_combined.txt'))
    e_grid_short /= e_grid_short @ pi_e
    n_e = len(pi_e)

1. **Their `ymarkov_combined.txt` is a generator `Q`, not a transition matrix** —
   the code applies `linalg.expm()` to it. Our `income_process_P.txt` is *already*
   `expm(Q)`. Dropping our `P` into that slot computes `expm(P)` — no error, just
   wrong numbers. **Fix:** export the combined generator instead. Because our
   combined chain is `P = kron(expm(Q₁), expm(Q₂)) = expm(kron(Q₁,I) + kron(I,Q₂))`,
   the right object to write is `Q_combined = kron(Q₁, I₂) + kron(I₁, Q₂)`, and
   their `expm` reconstructs our `P` exactly. Alternatively bypass their `expm`
   and load `P` directly. Either way, **add assertions on the matrix actually
   used**: rows sum to 1 and entries ≥ 0 (for `Pi_e`), or rows sum to 0 and
   off-diagonals ≥ 0 (for a generator).

2. **Mean normalisation is compatible.** They do `e_grid /= e_grid @ pi_e` —
   mean of `e` under the *stationary* distribution set to 1. Our `export_chain`
   does the same (`e / (pi @ e)`). Assert `pi_e @ e_grid ≈ 1` after load anyway,
   in case their row re-normalisation perturbs `pi_e`.

3. **State count is not hard-coded** — `n_e = len(pi_e)`, derived from the file.
   Our 75 states flow through; with `n_a = 200` the household grid is 75×200 vs
   KMV's 33×200. Slower, not breaking. (§14 keeps 33 available as a fallback.)

4. **Frequency is not stated in the code.** `r = 0.005` per period and "Income
   process from KMV 2018" both point to **quarterly** (KMV is quarterly; r≈0.5%/q
   ≈ 2%/yr). Our chain is quarterly (`P` = one-quarter transition), so if the
   model is quarterly we are aligned. **Verify explicitly in step 2** — check the
   discount-factor calibration target, the IRF horizon labels, and whether MPCs
   are reported quarterly or annual. If annual, feed `P⁴` (or re-discretize at
   annual and drop the §8d component-1 fix, since annual sampling carries the
   fast transitory component fine).

5. **`beta` is calibrated, not fixed** (`make_betas` + a separate calibration
   notebook, solved to an asset/wealth target given the income process). Swapping
   our chain **will change the calibrated `beta`**. The wealth-distribution
   comparison in the sanity path must therefore be done either at each country's
   own re-calibrated `beta`, or at a common fixed `beta` — a step-2 choice to
   make and document. This is expected, not a problem: it is the "reproduce
   steady state → swap → compare" plan working as intended.

6. Other calibration (one-asset HA block): `eis = 1` (log utility), `min_a = 0`
   (borrowing limit at zero), `max_a = 4000`, `n_a = 200`, `zeta = 0` (no
   cyclical income risk by default), `wN_aftertax = 0.7`. Nothing here assumes a
   particular income process beyond the chain and its mean-1 normalisation.

---

## 17. Task 5: is the iMPC matrix a sufficient statistic here? (checked against the source, not assumed)

**Question put by the supervisor.** Auclert–Rognlie–Straub ("The Intertemporal
Keynesian Cross", NBER WP 25020 / Econometrica) claim the matrix of
intertemporal MPCs (iMPCs) is a sufficient statistic for general-equilibrium
aggregate effects. If that holds for our model, the iMPC matrix — not the
direct/indirect decomposition — is the right object to compare across the
estimated income processes. Read the paper (full text extracted from the NBER
PDF via pypdf, since the standard PDF reader choked on it) rather than taking
the claim on faith.

**The exact theorem (Prop. 1, their eq. 13).** In their baseline "IKC
environment" (paper's Sections 2–5):
`dY = dG − M·dT + M·dY`, where `M` has entries `M_ts = ∂C_t/∂Z_s` (the
response of aggregate consumption at `t` to anticipated after-tax income at
`s`), and `M` alone is a sufficient statistic for the output path `dY` given
*any* bounded fiscal shock `{dG, dT}` satisfying the government's
intertemporal budget constraint (their eq. 12). Existence/uniqueness (Prop. 2)
requires `K(I−M)` invertible, `K` a known fixed operator.

**The conditions this requires — read directly off their Section 2.2–2.3 setup
(quoted/paraphrased from the extracted text):**
1. **A single asset** (real bonds only); no capital, no equity, no portfolio
   choice.
2. **Linear production, flexible prices, perfect competition** — real wage
   `w_t ≡ 1` always; sticky *wages* only. Output is purely demand-determined,
   `Y_t = C_t + G_t`, with no supply-side block to solve separately.
3. **Monetary policy pegs the real rate**: `1+i_t = (1+r)(1+π_{t+1})` forces
   `r_t = r` (constant) for *every* fiscal shock. Their own words: "the path
   `{r_s}` is not included as an argument [of the consumption function],
   because of our assumption that monetary policy fixes `r_s = r`." This is
   an **assumption about the policy closure, not a result** — it is the
   reason the theorem is stated only for *fiscal* shocks `(dG, dT)`, never for
   a shock that itself moves `r` (i.e., not for a monetary-policy shock).
4. Linearization around a steady state (Fréchet-differentiability of the
   consumption functional; a first-order/sequence-space-Jacobian result, not
   exact for arbitrary shock size — though the paper notes in a footnote this
   extends to first-order perturbations of a stochastic economy, not just
   deterministic MIT shocks).

**Extension to interest-rate/capital-gains shocks (Section 6, Props. 8/12) —
this is the part that matters for us.** They relax the fixed-`r` assumption
and add a second asset (equity) with valuation effects. The clean algebraic
result breaks, but a generalization survives: under **two extra conditions —
(a) EIS `σ = 1` (log utility) and (b) equal portfolio shares across agent
types** — the interest-rate Jacobian `M_r` is pinned down by `M` plus a single
extra object, `mcap` (the MPC out of unanticipated capital gains). Explicitly
covers their "HA-one" model (one-asset heterogeneous-agent with
discount-factor heterogeneity — i.e., our model family). They are explicit
that in this broader environment (their "quantitative environment", Section
7: capital, sticky prices, a Taylor rule) the sufficiency holds **numerically
in their calibration, not as a closed-form identity** like Prop. 1.

**Verdict for our model (ARS one-asset HANK, `hank_step2/`).**
- Our model uses a **Taylor rule reacting to inflation** (`kappa = 0.01`
  quarterly NKPC slope, §16) and runs **monetary-policy shocks** (an
  interest-rate shock *is* the experiment, STEP2_LOG §S1/S4/S5) — this is
  exactly the case excluded from Prop. 1 by construction (condition 3 above
  requires `r` fixed *regardless* of the shock, which is false for a monetary
  shock by definition). **Prop. 1's clean "`M` alone suffices" result does not
  apply to our headline exercise.**
- The **deficit-financed tax-cut** IRFs we already compute (S1 baseline) are
  closer in spirit to the paper's own object, but our model also carries a
  Taylor rule and a capitalization block there (not a real-rate peg), so even
  that exercise sits in their "quantitative environment", not their exact-
  theorem "IKC environment" — same caveat applies, with the same fallback
  below.
- **The generalization (Props. 8/12) plausibly does apply, conditionally.**
  Our model uses **log utility (`eis = 1`)** throughout (§16) and is
  **one-asset** with no cross-type portfolio heterogeneity (every household
  type holds the same combined bond/equity fund, so "equal portfolio shares"
  holds trivially) — both of the conditions Props. 8/12 require. So the
  *generalized* claim — that `M` plus `mcap` (not `M` alone) are sufficient
  for the interest-rate response, as established **numerically**, not as a
  closed-form theorem, in their own calibrated one-asset model — is the
  honest reading for us: **plausible, conditionally supported by the paper's
  own evidence for the same model class, but not a proven theorem for our
  specific calibration**, and explicitly not exact once a Phillips curve and
  Taylor rule are in play.
- **A separate, important distinction: sufficiency for *aggregate output* is
  not the same claim as sufficiency for the direct/indirect *decomposition*.**
  Prop. 1 is a statement about `dY` as a whole; it says nothing about how to
  split `dY` into channels. ARS's own Fig. 3b-style decomposition (which our
  §S3/S4 decomposition mirrors) is a *separate* accounting exercise built on
  partial-equilibrium counterfactuals (holding wages/transfers/capital-gains
  fixed one at a time) — it uses the same sequence-space-Jacobian machinery
  `M` comes from, but it is not what Prop. 1 is a theorem about.

**Bottom line.** The supervisor's intuition is directionally right but needs
one correction: it is not literally "the iMPC matrix `M`" that is the
sufficient statistic for our monetary-shock comparison — it is the **pair
`(M, mcap)`** (consumption response to income, and to capital gains), under
the log-utility + one-asset conditions our model satisfies, and the evidence
for that pairing being sufficient is numerical in ARS's own paper, not an
exact theorem the way Prop. 1 is for fixed-`r` fiscal shocks. With that
correction, the supervisor's suggested comparison is well-founded: compute
`M` (and `mcap`) for each estimated chain (KMV baseline, USA/FRA/GER/POL) from
the household block's Jacobian — SSJ produces this as an intermediate object
when solving steady state/Jacobians, so no new machinery is needed, only
extracting and tabulating it — and treat *that* as the primary cross-country
comparison object, with the direct/indirect decomposition kept as the
output-facing, more interpretable complement rather than replaced by it.

**Next step (deferred, not yet run as of this writing).** Extract `M`
(`∂C_t/∂Z_s`) and `mcap` from the household block's SSJ Jacobian for each of
the chains already swapped in step 2 (`hank_step2/swap_chain.py` machinery),
plot/tabulate the first several rows/columns (cf. their Figure 1) across
countries, and compare against the KMV/ARS baseline. This lives in
`hank_step2/`, a separate environment from this repo's `kmv_earnings/`
package.

---

## 18. Task 1: Poland estimated (2026-10-12)

**Data.** GRID Poland export (`data/dane grid PL.csv`), a *different* column
schema than the FRA/GER/USA `Stats_*.csv` (`std_perm`/`_resid_dNyr` instead of
`std_log_inc`/`_res_Nyr_log_chg`) — `grid_loader.py`'s `targets_from_grid_stats_csv`
now auto-detects both schemas (`_SCHEMAS`) and has a gender filter with an
assertion that exactly one gender survives (previously it silently averaged
across genders if a file had more than one — not a live bug for the
Male-only files used so far, but a latent one the task asked to close).
Male, prime age 25-55, **full native 2000-2020 window** (not the common
2001-2016 window — deliberately, see §0's open item above on the tension this
creates with the shape-comparison framing). Targets in `targets/pol_grid.json`.

**Estimation.** Same pipeline as the other three countries
(`scripts/close_step1.py`, now with `"POL"` added), full 6 free parameters
(no `beta2` pin — unlike GER/USA's canonical). DE and TikTak **agree** on the
persistent component (max rel. diff 1.3%, well under the 25% threshold used
for USA/FRA/GER) — no identification red flag from this check. DE is
canonical (OOS f=0.0037 vs TikTak's 0.0037, DE picked on a hair).

| param | value |
|---|---|
| lambda1 | 0.1623 |
| beta1 | 2.846 |
| sigma1 | 1.979 |
| lambda2 | 0.00877 |
| beta2 | 0.0169 |
| sigma2 | 1.401 |

Implied ergodic sd of the persistent component: **0.71**, below KMV's 0.95 and
below every other country here (USA free 1.22, FRA 1.09, GER free ~2.68) —
Poland's persistent component is the *most* mean-reverting of the four, not
the least, on the free/unrestricted estimate. `beta1 = 2.85` sits closer to
USA (2.75-2.97 depending on restriction) than to FRA/GER (4.3-5.5) — on
transitory persistence Poland looks more USA-shaped, which is the dimension
DECISIONS' two-asset trigger (bottom of the open-decisions list) cares about;
on persistent-component dispersion it looks nothing like the GRID-USA problem
case (§12). No parameter at a bound. Fit is good across all 8 moments
(largest relative error 4.7%, on `frac_d1_lt_10`). Discretized chain shows the
same known gaps as the other countries (var_d1 ~27% low, §8d; kurtosis
overshoots a bit at this grid size) — not new issues, no new mechanism.

**Schema-mapping verification (2026-10-12, requested alongside the two-asset
trigger check below).** The Polish GRID file uses different column names than
the FRA/GER/USA file (`std_perm`/`_resid_dNyr` vs `std_log_inc`/`_res_Nyr_log_chg`),
handled by `_SCHEMAS` auto-detection in `grid_loader.py` (§ above). Checked that
this is a *semantically* correct mapping, not just one that parses: for 2010,
computed `var_log_earns`, `var_d1`, `kurt_d1`, `var_d5`, `kurt_d5`, and all
three change-fractions **by hand** directly from the raw CSV row, hard-coding
the percentile-probability grid myself rather than reusing any of
`grid_loader.py`'s column-detection code — a genuinely independent path, not
just a second call into the same helper functions. Matched the loader's
output exactly on every quantity (e.g. `frac_d1_lt_10` hand-computed 0.4919 vs
loader 0.4919). Also confirmed `std_resid_d1yr` and `std_resid_d5yr` (and
their kurtoses) are genuinely different columns feeding the right targets, not
an accidental same-column reuse. No mapping error found.

**For the write-up.** Poland does not resemble the GRID-USA failure case from
§12/STEP2_LOG §S5 (thin tail, slow-decaying transitory component) — its
`kurt_d1` (15.24, data) and `beta1` sit in between the FRA/GER and USA
patterns. Whether this is close enough to FRA/GER-shaped for the one-asset +
death HANK structure to calibrate is a step-2 question, not answered here —
see the two-asset trigger note.

---

## 19. Task 3: bootstrap standard errors (2026-10-12)

**Design.** `scripts/bootstrap_se.py`. Resample years with replacement from
each country's own GRID year panel (within its estimation window), recompute
the 8 target moments as the simple across-year mean (verified bit-for-bit
identical to the original `targets/<cc>_grid.json` construction when applied
to the unresampled window), re-estimate with a warm start at the country's
point estimate. 200 reps/country (800 total), DE only (not TikTak — §9 already
establishes TikTak's local phase is unreliable here, and it would roughly
double the cost for no identification benefit). GER/USA keep their `beta2`
pin in every rep, so their reported SE/CI cover only the five free
parameters — bootstrapping a calibrated quantity would misrepresent
calibration as sampling uncertainty.

**Budget.** Full-serious settings (50k-worker panel, DE maxiter=100/popsize=12,
NM polish 400 iters) cost ~30 min per estimate even at 10-way parallelism —
200x4 reps at that cost would run for weeks. Reduced to a 5,000-worker
panel, DE maxiter=10/popsize=5, NM polish 50 iters, still warm-started at the
point estimate. Validated against a larger reference budget (15,000 workers,
maxiter=20/popsize=6, NM 100 iters) on 3 resamples per country before
committing to the full sweep, as the task required — not skipped.

**A methodological incident, worth recording because it could recur.** The
first validation run flagged what looked like a severe reduced-budget bias
specific to Poland (`beta1` off by 40-67%, nothing like USA/FRA/GER's <10%
gaps). Root cause: the validation step launched **while Poland's serious
`close_step1.py` run was still mid-flight**, so it warm-started from a stale
`output/pol/params.json` left over from an earlier `--quick` smoke test
(which, carelessly, wrote to the same path as the serious run) — not from the
real estimate. Confirmed by checking the file's mtime and content against the
quick-test log. The real 800-rep sweep was unaffected in the end, because its
jobs are queued country-by-country (Poland last) and by the time Poland's
chunk ran, the real file was in place — but the validation's **aggregate**
verdict ("OK, not systematically biased") was also wrong on its own terms: it
pools the systematic-bias check across all four countries, so Poland's large
deviation was averaged away by USA/FRA/GER's small ones (pooled
`beta1` mean=-11%, sd=27%, under the 15%-mean threshold; the per-country
numbers were -39% to -67% for Poland alone). Fixed `report_validation` to
check per-country, not pooled, so one bad country can no longer hide behind
the others. **Lesson for later reduced-budget/warm-start exercises in this
project (Task 4 included): never launch a bootstrap/recovery/reduced-budget
job against a point estimate that is concurrently being overwritten by
another running job, and never write quick-smoke-test output to the same path
as a serious run's output.**

**Result — final 800-rep sweep (post-fix, Poland using the correct warm
start), SE and 95% percentile CI:**

| | lambda1 | beta1 | sigma1 | lambda2 | beta2 | sigma2 |
|---|---|---|---|---|---|---|
| USA mean | 0.214 | 2.968 | 1.942 | 0.01045 | 0.01398 [locked] | 1.561 |
| USA SE | 0.0040 | 0.082 | 0.0076 | 0.00031 | - | 0.021 |
| FRA mean | 0.190 | 4.367 | 1.915 | 0.00694 | 0.00484 | 1.094 |
| FRA SE | 0.0022 | 0.141 | 0.0096 | 0.00023 | 0.00040 | 0.0029 |
| GER mean | 0.156 | 5.518 | 1.798 | 0.01029 | 0.00663 [locked] | 1.083 |
| GER SE | 0.0018 | 0.126 | 0.0067 | 0.00009 | - | 0.0021 |
| POL mean | 0.159 | 2.850 | 1.982 | 0.00882 | 0.01690 | 1.411 |
| POL SE | 0.0033 | 0.055 | 0.0113 | 0.00035 | 0.00119 | 0.0160 |

Poland's bootstrap mean (`beta1`=2.850, 95% CI [2.74, 2.96]) sits right on the
point estimate (2.846) with a tight CI — direct confirmation the earlier
validation scare was the warm-start artifact above, not a real problem with
the reduced budget. All four countries' point estimates sit comfortably
inside their own bootstrap means, and relative SEs are small (mostly 1-5% of
the point estimate) except where a parameter is already known to be weakly
identified over a short panel.

**Honest caveat (as the task required, stated plainly, not buried).** This
measures year-to-year *sampling* variability in the reported GRID statistics
only. It is not a full-sample-uncertainty SE — GRID's underlying microdata and
true sampling distribution are not available, so this cannot capture, e.g.,
measurement error in the underlying earnings records or the finite-population
correction a full survey-design SE would need.

**Second caveat (added 2026-10-12 after §21): the CIs are conditional on the
correct basin.** Every bootstrap rep is warm-started at the country's point
estimate, so by construction each rep searches locally around the basin the
point estimate already sits in. That is the right design for measuring
sampling variability cheaply, but it means the reported SE/CI contain no
component for the risk that the point estimate itself sits in the wrong
basin. §21 shows that risk is not hypothetical: for GER, a blind
re-estimation from a neutral start reaches a different point in 4 of 4
draws, and for POL in 1 of 4. So GER's `β₁` SE of 0.126 should be read as
"variability of the estimate given this basin", not as total uncertainty
about `β₁`. Say this explicitly in the thesis wherever the SE table appears.

**For the write-up.** Report the SE/CI table with the caveat above stated in
the same sentence, not in a footnote. The validation incident is worth one
sentence as a methods point (warm-started reduced-budget resampling needs a
*stable* warm start) but not a deep dive.

---

## 20. Poland's two-asset trigger check, and whether it depends on the window choice (2026-10-12)

**Question 1: does Poland need two-asset, or does one-asset + death suffice?**
Run `hank_step2/swap_death.py pol` (one-asset ARS model + ARS β-heterogeneity
+ KMV stochastic death, jointly re-calibrated to `A=20`, impact labor
`MPC=0.20`, SCF Lorenz — identical procedure to FRA/GER/USA, STEP2_LOG §S5).
**Poland calibrates cleanly**, matching FRA/GER, not the GRID-USA failure
mode:

| | target | **POL** | FRA | GER | USA (degenerate) |
|---|---|---|---|---|---|
| converged | — | **yes**, resid ~1e-12 | yes | yes | no |
| MPC (labor, impact) | 0.20 | **0.2000** | 0.200 | 0.200 | 0.024-0.029 |
| frac at borrowing constraint | — | **0.226** | 0.225 | 0.250 | ~0.01 |
| dβ (heterogeneity spread) | — | **0.125** | 0.093 | - | →0 |
| β_lo | — | 0.873 | 0.905 | - | - |
| decomp yr1 direct/indirect | — | **22.6/77.4%** | 25.5/74.5 | 24/76 | degenerate |
| monetary dY impact (base 2.41) | — | 2.30 (-4.5%) | 2.29 (-5%) | 2.33 (-4%) | 1.81-1.85 |

Every diagnostic that separates the FRA/GER cluster from the USA failure
(dβ staying well off zero, a believable hand-to-mouth fraction, "indirect
dominates" surviving intact) lands on the FRA/GER side for Poland. **No
escalation to two-asset is forced by this test.** Note this contradicts the
rule-of-thumb threshold recorded in the open-decisions list (`β₁ ≲ 3` →
GRID-USA-shaped): Poland's `β₁ = 2.85` sits just below it, but the ergodic
variance of its persistent component (sd 0.71, the *lowest* of all four
countries — §18) and the shape of its kurtosis pull the other way, and the
actual calibration run is unambiguous. **The lesson: run the check, don't
apply the threshold mechanically** — exactly the caution raised before running
it (DECISIONS §0's two-asset trigger note).

**Question 2: does this answer depend on the 2000-2020 vs a split window?**
Motivated by the 2009 kurtosis/small-change-fraction regime shift documented
in `output/pl_diagnostics/README.md` (kurt Δ1y jumps from ~13.5 to 15-17
starting exactly in 2009 and stays elevated for seven years, at roughly
unchanged variance). Built separate targets for 2000-2008 and 2009-2020
(`targets/pol_grid_{pre,post}2009.json`, `scripts/pol_regime_split.py`,
DE-only — TikTak not needed for this exploratory split, §9), then ran the
same `swap_death.py` check on both resulting chains.

*Parameters move only modestly across the split* (full-window shown for
reference):

| | full (00-20) | pre (00-08) | post (09-20) |
|---|---|---|---|
| lambda1 | 0.162 | 0.177 | 0.153 |
| beta1 | 2.846 | 2.857 | 2.816 |
| sigma1 | 1.979 | 1.986 | 1.976 |
| lambda2 | 0.0088 | 0.0096 | 0.0078 |
| beta2 | 0.0169 | 0.0178 | 0.0130 |
| sigma2 | 1.401 | 1.397 | 1.333 |
| ergodic sd(z2) | 0.714 | 0.726 | 0.729 |

Despite the sharp raw-moment shift (kurt_d1 14.3 -> 16.0, var_d1 0.285 ->
0.239 across the split — moving in *opposite* directions, i.e. a genuine
shape change, not a scale change), the estimated process absorbs it with a
fairly small change in the persistent component's activity (`lambda2`/`beta2`
both a bit lower post-2009 — slightly rarer, slightly slower-decaying
persistent shocks), not a different process family. Both sub-windows fit well
(f=0.0046 pre, f=0.0032 post, vs f=0.0037 full-window).

*And the two-asset trigger answer is identical on both sub-windows*:

| | full | pre-2009 | post-2009 |
|---|---|---|---|
| converged | yes | yes | yes |
| MPC | 0.200 | 0.200 | 0.200 |
| frac HtM | 0.226 | 0.220 | 0.225 |
| dβ | 0.125 | 0.141 | 0.112 |
| decomp yr1 direct/indirect | 22.6/77.4 | 22.6/77.4 | 22.9/77.1 |

**Answer: no, the two-asset decision does not depend on the window choice.**
All three windows calibrate cleanly, with near-identical MPC, hand-to-mouth
share, and decomposition. This is itself worth reporting as a robustness
result, not just a technicality dispatched in a footnote: the 2009 regime
shift is real and visible in the raw moments (§ diagnostics), but it is not
large enough, nor of the right *kind*, to move the one-asset-vs-two-asset
decision. Had the answer differed between sub-windows, the choice of
2000-2020 vs a split would have become a substantive modelling decision
rather than a data-availability one (§0's open item) — it did not.

**Schema-verification note (requested alongside this check, done, see §18).**

**For the write-up.** This section answers the two-asset trigger question
definitively for Poland and closes the corresponding open-decisions bullet
above. Report the full-window numbers as primary and the pre/post split as a
robustness check demonstrating the window choice is immaterial to this
decision.

---

## 21. Task 4: parametric recovery — a sharper identification test than the line scan (2026-10-12)

**Design (reduced scope, see the time-estimate correction below).**
`scripts/recovery_exercise.py`. Take each country's FREE (unrestricted)
6-parameter estimate as the "true" process — for GER this required a new
free estimate first (§ above, `output/ger_free_beta2/`), since GER's canonical
pins `beta2`. Simulate a synthetic panel from the true process with a
generating seed **different** from the CRN seed used inside the estimator
(otherwise the test is trivial — the "recovered" estimate would just read
back the same random draws it was fitted to). Compute the 8 moments on that
panel, treat them as data, and re-estimate from scratch with DE (not
TikTak — §9), warm-started at the *neutral* KMV Table-3 values, not at the
truth (warm-starting at the truth would make even a weakly-identified
parameter trivially "recover" and defeat the point). Compare recovered vs
true parameters, especially `beta2` and the implied ergodic sd of the
persistent component, across several generating seeds and panel sizes.

Reduced to **2 seeds x 2 panel sizes (50k, 300k)** = 16 reps total (not the
original 3x3=36) after the first timed rep showed ~28-30 min per
re-estimation even at full 10-core parallelism (the initial ~10 min guess,
based on under-calibrated assumptions about DE-only cost, was wrong by
~3x) — would have cost ~17-18h at full scope; the reduction was agreed with
the user to fit an ~8h budget. Also ran concurrently with the Task 1/2
two-asset trigger checks (§20) and Task 2's nested-specs runs, at the user's
explicit instruction to accept the resulting slowdown.

**Headline result, properly disentangled (optimiser reliability vs
identification — the same distinction §9 insists on, now with sharper
evidence):**

| country | beta2 reldiff (mean, 4 reps) | consistent across reps? | improves with panel size? |
|---|---|---|---|
| USA | -1.8% (range -3.1% to -0.6%) | yes, tight (sd 1.0%) | no (if anything slightly worse at n=300k: -0.9%->-2.7%) |
| FRA | -4.0% (range -4.4% to -3.5%) | yes, tight (sd 0.3%) | no (-3.7%->-4.3%) |
| **GER** | **-82.3%** (range -80.7% to -83.2%) | **yes, every single rep** (sd 1.0%) | **no** (-81.9% at both panel sizes) |
| POL | -25.2% mean, but **3/4 reps within 1%**, 1/4 catastrophic (-99.4%) | **no — bimodal** | the failure was at the *larger* panel (n=300k), not smaller |

**USA and FRA: beta2 is well identified.** Small, stable biases (1-4%), no
dependence on sample size — ordinary small-sample bias, not an
identification problem.

**GER: beta2 (and, surprisingly, beta1/lambda1 too) fails systematically —
real evidence, but it is evidence of a flat valley, not two equally good
optima.** All 4 independent generating draws converge to the **same**
alternate point (`beta1~1.84, lambda1~0.09, beta2~0.0001`) regardless of seed
or panel size — the consistency itself (not scatter) is the signal: this is a
robust, repeatable feature of the objective surface for GER-shaped data, not
optimiser noise. But checking the objective value at the true point vs the
recovered point on the *same* synthetic targets shows the true point is
genuinely better (`f=0.00024` at truth vs `f=0.00085` at the recovered point,
on one rep) — both are small, "good fit" numbers by this project's own
standard, so the surface is very flat between them, but it is not literally
indistinguishable. **Reading: with a fully serious optimisation budget — the
exact DE configuration used for every real-country estimate in this
project — the search still cannot reliably climb out of a shallow secondary
basin back to the true point, for a process this weakly identified in its
persistent component.** That is a stronger and more operational statement
than "the objective sweep is flat near the optimum" (§8b): it says the actual
*estimation procedure*, not just the *objective function in the abstract*,
fails to recover the true process here.

**A striking cross-check: this alternate basin is not new.** `beta1 ~ 1.84` is
in the same low-beta1 neighbourhood TikTak's Nelder-Mead local phase got stuck
in on GER's *real* data (§9: `beta1 ~ 2.2-2.4`) — a completely different
optimiser (sequential Nelder-Mead vs differential evolution), a completely
different dataset (real GER moments vs noise-free synthetic data generated
from the model's own free estimate), and they land in the same region. That
is strong independent evidence this is a genuine, repeatable feature of the
objective landscape for GER/near-random-walk-shaped data, not an artefact of
one algorithm or one dataset.

**Does this put the canonical free-GER estimate itself in doubt?** Checked
directly: evaluated the "alternate basin" point found in the recovery
exercise against GER's **real** (not synthetic) targets. `f=0.00564` at the
alternate point vs `f=0.00219` at the canonical free estimate — the canonical
estimate is still clearly better (2.6x) on the real data, a wider margin than
on the synthetic data (3.5x on one rep, but both numbers tiny). **So the free
GER estimate used elsewhere in this project is not undermined** — but the
margin is not large, and this is exactly the kind of check that should be
repeated if GER's free estimate is ever leaned on harder (e.g. if the pinned
`beta2` calibration, §8b, were ever revisited).

**Poland: mostly well identified, but not bulletproof — a genuine,
occasional optimiser failure, distinct from GER's systematic one.** 3 of 4
reps recover `beta2` to within 1%. The 4th (seed 202, n=300,000) collapses to
the *same qualitative* alternate basin as GER (`beta2 -> 0.0001`,
`lambda2` exploding +1039%, ergodic sd blowing up to 4.7 vs true 0.71). But
checking the objective at the true point vs the recovered point on that
rep's synthetic targets shows a **large** gap (`f=0.00011` at truth vs
`f=0.00947` at the recovered point — 83x worse, not 3.5x as for GER). That
asymmetry is the tell: **this is DE failing to converge on one unlucky run,
not a robust alternate optimum** — if it were a genuine competing basin like
GER's, it would show up with a comparably small objective gap, and ideally in
more than 1 of 4 draws. Still informative: it shows the same qualitative
trap exists within reach of Poland's parameter region too (never triggered
for USA or FRA in 8 combined reps), so Poland's identification margin is
real but **not as comfortable as USA/FRA's** — multiple-seed or multi-start
confirmation (as already done for the canonical estimate via DE-vs-TikTak
agreement, §18, which agreed to 1.3%) is worth keeping as standard practice
for Poland specifically, not dropped as redundant.

**For the write-up.** This is the headline identification result for the
thesis, stronger than the line scan (§9) because it does not require already
knowing where to look — it is a blind re-estimation from a neutral start.
State the three-way distinction plainly: USA/FRA well-identified;
**GER systematically fails even with a full serious budget and a
noise-free, best-case synthetic dataset** — the clearest, most direct
evidence in the whole project that `beta2` cannot be estimated for a
near-random-walk persistent component from these 8 moments, which is exactly
why it is calibrated rather than estimated (§8b); **Poland is well identified
on average but has a real, if infrequent, failure mode** worth flagging as a
caveat on the point estimate's robustness, not as grounds for recalibrating.

---

## 22. Task 2: nested specifications — what does the jump-drift structure add? (2026-10-12)

**Design.** Three restrictions nested inside the full 6-parameter process,
estimated for all four countries (`scripts/nested_specs.py`, DE-only, serious
budget), compared on the full 8-moment fit against each country's own free
6-parameter estimate:
- **V1, single component** (3 free params): component 2 switched off
  (`sigma2=0`), Poisson arrival kept for the surviving component.
- **V2, deterministic arrival** (2 free params): V1 plus shocks arriving
  every sub-period with probability 1 instead of at a Poisson rate — KMV's
  own footnote nests a discrete-time AR(1) in logs as this limit.
- **V3, `β₂ = 0`** (5 free params): persistent component made a pure random
  walk (no mean reversion), everything else free.

**V1 (single component): destroys the LEVEL variance specifically, leaves
shape and 1-year dynamics almost untouched — not what a naive guess would
predict.** `var_log_earns` collapses to roughly a third to a half of its data
value for every country (USA 0.956→0.363, FRA 0.476→0.193, GER 0.645→0.163,
POL 0.508→0.284) — the single biggest driver of V1's 85-125x worse objective.
But `var_d1` is barely touched (within 5-7% of data, all four countries), and
`kurt_d1` is only mildly worse (USA 12.8→11.8, FRA 15.2→15.2 essentially
unchanged, GER 17.8→18.1, POL 15.2→14.4). `kurt_d5` and `var_d5` degrade
moderately (the second component's own occasional jumps matter at the
5-year horizon).

**Reading (corrected 2026-10-12): one component cannot simultaneously fit
level dispersion and short-run change dynamics.** The original framing here
("the second component's job is not primarily to add leptokurtosis, its job
is to carry level dispersion") overclaimed a structural division of labour
that the evidence does not establish. With only 3 free parameters chasing 8
equally-weighted moments, *which* moment gets sacrificed is an artefact of
the SMM objective's equal weighting, not proof of what a second component is
structurally "for". A differently-weighted objective (e.g. upweighting
`var_log_earns`) could plausibly push V1 to sacrifice `kurt_d1` or `var_d5`
instead and preserve the level variance better. State the finding as *what
breaks under this particular objective*, not as a structural division of
labour between the components.

**V2 (deterministic arrival): uniformly catastrophic, and the kurtosis
collapse is exact — a clean confirmation of KMV's footnote.** `kurt_d1` falls
to **3.00-3.04 for all four countries**, the Gaussian value to three
significant figures. **This part is mechanical, not a finding** — summing
Gaussian innovations every sub-period, with no Poisson gaps between them, is
itself a Gaussian process by construction (a sum of Gaussians is Gaussian),
so kurtosis landing at exactly 3 is guaranteed by the restriction's
definition and confirms nothing beyond "the arithmetic was implemented
correctly". **The informative part is everything else**: every other moment
is also far off (`frac_d1_lt_10` roughly a third of its data value
everywhere, `var_log_earns` similarly depressed) — i.e. deterministic arrival
does not merely erase kurtosis while leaving the rest of the distribution
intact, it degrades the *entire* moment profile. That breadth, not the
(guaranteed) kurtosis collapse, is the evidence that genuine Poisson timing
is doing real work.

**V3 (`β₂ = 0`, random walk): costs little in FIT for every country, but
"little cost" hides a very different story for GER than for USA/FRA/POL.**
OOS objective: USA 0.0058→0.0212 (3.6x), FRA 0.0013→0.0108 (8.3x), POL
0.0037→0.0130 (3.5x), **GER 0.0025→0.0028 (1.1x — and per §9's Monte Carlo
noise check, a 1.1x ratio at this objective scale is within the ~23%
simulation-noise band, i.e. indistinguishable from "no cost at all")**.
Every individual moment stays close to the free estimate for every country —
`var_log_earns` within 1-3%, `kurt_d1` within 2-4%, fractions within a few
percentage points.

**Correction (2026-10-12): "every moment close" hides a complete change of
process for Poland, and the 3.5x penalty understates how costly V3 really
is for PL.** Poland's V3 estimate sits at or against its parameter bounds:
`λ₂ = 0.09999636` (bound `[0.001, 0.10]` — at the ceiling), `σ₂ = 0.2129`
(bound `[0.20, 3.5]` — against the floor), and `β₁ = 0.50` where the free
estimate has `β₁ = 2.85` — a more than 5-fold change in the *transitory*
component, not just the persistent one. The optimiser is not finding "the
same process with `β₂` turned off"; it is finding a **structurally
different process** (a much slower-decaying transitory component, pushed up
against two bounds on the persistent side) that happens to reproduce the 8
moments about as well. An objective ratio that looks mild (3.5x) is
therefore not a fair measure of how much this restriction costs Poland
economically — the parameter vector itself has moved into a different
region of the space entirely. USA and FRA's V3 estimates were not checked
against their bounds here and should be before repeating this claim for
them.

**Synthesis with Task 4 (§21) — corrected (2026-10-12): the original
synthesis here overstated the generalization and directly contradicted
§21.** The original text claimed `β₂` is "only weakly disciplined... in
general," with GER "merely the most extreme case of a pattern that holds...
for all four countries." That is not what §21 shows: USA and FRA recovered
`β₂` from synthetic data to within **1-4%** across every draw, and Poland
recovered it to within **1%** in 3 of 4 draws. Those are well-identified
results, not weak ones. The correct, narrower synthesis: **GER's `β₂` is
genuinely unidentified** (§21: systematic recovery failure in 4/4 draws,
with a known alternate basin; this section: a 1.1x fit cost that is itself
noise-level). **USA, FRA, and Poland's `β₂` are identified** — §21 recovers
them well — **but fixing `β₂` to zero still costs only a small absolute
amount of fit for them too** (3.5-8.3x on an `f` that starts below 0.006),
*and*, per the correction above, the parameter vector that achieves that
small cost for Poland is economically a different process, not a
free lunch. **Whether this small-in-objective-terms restriction matters for
anything that is actually used downstream (the HANK calibration) is a
separate question from whether it is "identified" in the SMM sense — §23
is the test that actually answers it**, by profiling the ergodic variance
`β₂` maps into rather than treating the fit cost as decisive on its own.

**For the write-up.** Report the full 8-moment table
(`output/nested_specs/nested_specs_comparison.csv`) as a figure or table —
it is self-explanatory. The headline sentence, corrected: *the jump
structure is essential (V2 fails completely on every moment, not just
kurtosis — confirming KMV's AR(1) footnote); one component cannot fit level
dispersion and short-run dynamics simultaneously, though which one breaks is
an artefact of equal SMM weighting (V1); and the persistent component's
decay rate is weak-to-well identified depending on the country — genuinely
unidentified for Germany, identified but cheap-to-restrict for USA/FRA/
Poland, with Poland's "cheap" restriction masking a different process
entirely (V3, see §23 for whether any of this matters for HANK outcomes).*

---

## 23. Does beta2's weak fit-discipline actually matter for HANK outcomes? A Poland ergodic-variance profile (2026-10-12)

**Motivation.** §21/§22 show `β₂` is weakly disciplined by the 8-moment SMM
objective (unidentified for GER, cheap-to-restrict for USA/FRA/POL). But
`β₂` itself never enters the HANK model — what enters, via the discretized
chain's ergodic distribution and `death.py`'s newborn draws, is the
persistent component's **ergodic variance**, `var2 = λ₂σ₂²/(2β₂)`. The
question §22 could not answer on its own: does moving `var2` around, within
the range the SMM objective is roughly indifferent to, actually change
anything a HANK run would report? This is the test `scripts/pol_beta2_profile.py`
and the follow-up `swap_death.py`/`impc_matrices.py` runs were built to
answer directly, rather than inferring it from the fit cost alone.

**Design.** Pin `var2` via `ErgodicVarConstraint`, re-estimate the other 5
parameters at Poland's full canonical budget (50k-worker panel, DE
maxiter=100/popsize=12, NM 400, warm-started at `output/pol/params.json`,
`workers=10`), score OOS on the same 100k-worker/seed=987654 panel
`nested_specs.py` uses, export the chain, then run `swap_death.py` (one-asset
+ β-heterogeneity + KMV death, jointly calibrated to `A=20`, `MPC=0.20`, SCF
Lorenz — identical procedure to every other country in this project) and
extract the iMPC sequence. Grid: `var2 ∈ {0.30, 0.40, 0.70, 0.9103}`; Poland's
own canonical (free) estimate sits at `var2 = 0.5094` (ergodic sd 0.714).
`0.9103` is KMV's own anchor — the same ergodic-variance target already
pinned for USA and GER (§8b, and the USA pin in the open-decisions list) —
chosen deliberately as the row that determines whether a PL-vs-USA
comparison built on a common ergodic-variance anchor is clean.

A reduced-budget version of this exercise, run in a separate investigation
before this one, found `var2=0.30` breaks `var_log_earns` (model 0.375 vs
data 0.508) and `var2=0.9103` breaks `var_Δ5y` (model 0.410 vs data 0.477),
with every other moment close to data in both cases. **The full-budget run
below reproduces those numbers almost exactly** (OOS `f`: `0.0753` vs the
preliminary `0.075` at `var2=0.30`; `0.0310` vs `0.031` at `var2=0.9103`) —
the reduced budget was already adequate for this particular question, and
the full run adds the two middle points plus the HANK-side outcomes.

**The full table:**

| var2 | ergodic sd | OOS f | worst moments (rel. err) | converged | MPC | frac HtM | dβ | β_lo | decomp yr1 dir/indir | monetary dY impact | deficit dY impact |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.30 | 0.548 | 0.0753 | var_log_earns −26.3%, var_d1 +6.1%, frac_d1_lt_20 −3.8% | **NO** | 0.319 | 0.329 | **0.300 (bound)** | 0.700 | 19.4/79.9 | 2.546 | 0.524 |
| 0.40 | 0.632 | 0.0203 | var_log_earns −13.0%, var_d5 +3.2%, frac_d1_lt_10 +3.1% | yes | 0.200 | 0.253 | 0.120 | 0.879 | 21.0/79.0 | 2.279 | 0.419 |
| **0.509 (canonical, free)** | **0.714** | **0.0037** | (fit is good on all 8 moments, §18) | **yes** | **0.200** | **0.226** | **0.125** | **0.873** | **22.6/77.4** | **2.303** | **0.400** |
| 0.70 | 0.837 | 0.0200 | var_d5 −8.9%, var_log_earns +6.7%, frac_d1_lt_10 +6.2% | yes | 0.200 | 0.181 | 0.107 | 0.893 | 19.6/80.5 | 2.387 | 0.480 |
| 0.9103 (KMV anchor) | 0.954 | 0.0310 | var_d5 −13.5%, var_d1 +6.5%, var_log_earns +6.4% | **NO** | 0.333 | 0.340 | **0.300 (bound)** | 0.700 | 18.7/80.8 | 2.533 | 0.540 |

**Headline result: the calibration fails at BOTH ends of the tested range,
not just the low end — and KMV's own anchor is one of the failures.** At
`var2 = 0.30` and `var2 = 0.9103`, `β`-heterogeneity (`dβ`) is driven to its
upper bound (0.30) and `β_hi`/`β_lo` to the edge of their box, and the joint
`(A, MPC, Lorenz)` calibration still **overshoots** the MPC target
(0.319-0.333 instead of 0.200) rather than converging — a genuine
infeasibility, not a slow convergence. Only `var2` roughly in `[0.4, 0.7]` —
bracketing Poland's own data-implied value of `0.509` — supports a clean
calibration. The canonical free estimate is also, by a wide margin, the
best-fitting point in SMM terms (`f = 0.0037` vs `0.020-0.075` for every
pinned alternative), consistent with §22: pinning `var2` away from its own
optimum costs more than §22's bare `β₂=0` comparison suggested, once the
*direction* of the restriction (too low vs too high) is allowed to vary
rather than only checking `β₂=0`.

**Direct answer to the motivating question: yes, it matters, and it matters
specifically for the PL-vs-USA comparison.** `var2 = 0.9103` is exactly the
ergodic-variance anchor imposed on USA and GER (§8b) for exactly the reason
that their own free estimates do not support a clean calibration. Poland is
different: its own free estimate (`var2 = 0.509`) **does** support a clean
calibration — but forcing it onto the *same* anchor used for USA breaks that.
So a PL-vs-USA comparison conducted by pinning both to a common
ergodic-variance target is not available as a clean design for Poland the
way it is already accepted as necessary for USA — Poland's clean calibration
depends on using its own data-implied dispersion, not an imposed common one.
**This means part of any PL-vs-USA difference in the main comparison
reflects Poland running on its own data-implied dispersion while USA runs on
an imposed one, and that has to be stated in the main text, not left as an
appendix caveat** — exactly the risk flagged when this exercise was
commissioned.

**What moves across the profile, substantively.** `frac_htm` falls
monotonically from `0.33` (var2=0.30, not converged) through `0.25`, `0.23`
(canonical), `0.18` (var2=0.70) back up to `0.34` (var2=0.9103, not
converged) — U-shaped, mirroring the OOS fit. The monetary-shock impact
response and the direct/indirect decomposition stay qualitatively similar
across all five rows (direct share 18.7-22.6%, "indirect dominates" survives
everywhere) — the decomposition's *qualitative* conclusion is robust to
where `var2` sits, even though the *calibration* is not. iMPC sequences
(`hank_step2/results/impc_profile.png`) show the two non-converged points
visibly overshoot the impact response (0.32-0.33 vs 0.20) while tracking
each other closely beyond the impact quarter — consistent with the
overshoot being concentrated in the contemporaneous MPC, not the whole
sequence shape.

**For the write-up.** This is the section that should accompany any
PL-vs-USA HANK comparison, not §22 alone. State plainly: Poland's one-asset
+ death calibration is clean at its own data-implied ergodic dispersion and
breaks (fails to converge) if forced onto USA's imposed anchor; this is
itself a finding about how the countries differ, not a nuisance to be
normalized away.

### 23a. The OOS-f-vs-var2 plot, and which points are noise-indistinguishable

`scripts/plot_pol_var2_profile.py`, `output/pl_diagnostics/11_var2_profile.png`.
Plots OOS `f` (log scale) against `var2` for the 4 pinned profile points
plus the canonical free estimate, and marks a noise-indistinguishability
threshold at **+20% of the minimum `f`**, per §9's direct measurement that
Monte Carlo noise at this objective scale runs ~23% of the mean.

**Correction (2026-10-12): the threshold must be measured from the true
minimum — the canonical free estimate (`f=0.0037`) — not from the best
*pinned* point (`var2=0.70`, `f=0.0200`).** The first version of this plot
used the latter, which is the wrong reference: the pinned grid are
constrained re-estimates of a *different* (5-parameter) problem, and
comparing them only to each other cannot say anything about distance from
the actual optimum. With the correct threshold
(`1.20 × 0.0037 = 0.0044`), **none of the four tested pinned values falls
within noise of the true minimum** — every one of them, including the
"best" pinned point at `var2=0.70`, is a real, non-noise degradation
relative to free estimation. A dedicated **control point**, pinning `var2`
at exactly the canonical's own value (`0.5094`) and re-running the full
profile pipeline, is the direct test of whether the pin itself (as opposed
to pinning it at the *wrong* value) costs anything — it should reproduce
the canonical free estimate closely if the constraint is genuinely
non-binding there.

**Control point, run and confirmed (2026-10-12):** `OOS f = 0.00366` against
the canonical free estimate's `0.00367` — indistinguishable — and the
recovered parameters match the free estimate to 3-4 significant figures
(`λ₁=0.16233` vs `0.16232`, `β₁=2.8455` vs `2.8462`, `σ₂=1.4008` vs
`1.4010`, etc.). This validates both the profile pipeline itself (pinning a
non-binding constraint reproduces the unconstrained optimum, as it should)
and the choice of `f=0.0037` as the reference minimum for the noise band
above. Now included in `output/pl_diagnostics/11_var2_profile.png`.

**Among the pinned points, `var2=0.40` and `var2=0.70` remain
noise-indistinguishable *from each other*** (`0.0203` vs `0.0200`, well
within the ~23% noise band) — a real, narrower claim than the one above:
the data cannot tell these two *pinned* values apart from one another, but
both are nonetheless clearly worse than letting `var2` float freely.
`var2=0.30` and `var2=0.9103` are not close to either the `[0.4,0.7]`
pinned pair or the true minimum — real degradation on both sides, not
noise, under either reference point.

**Which moment breaks on each side — confirmed cleanly, and consistently on
both sides of the band, not just at the extremes.** `var2=0.30`'s worst
moment is `var_log_earns` (`-26.3%`); `var2=0.40`'s worst moment is *also*
`var_log_earns` (`-13.0%`, same direction, milder) — the low side is
disciplined by **level dispersion**, consistently. `var2=0.9103`'s worst
moment is `var_Δ5y` (`-13.5%`); `var2=0.70`'s worst moment is *also*
`var_Δ5y` (`-8.9%`, same direction, milder) — the high side is disciplined
by **5-year change variance**, consistently. This is exactly the asymmetric
discipline hypothesized before running the full-budget version, now
confirmed moment-by-moment on four independent points, not just the two
extremes: `var_log_earns` pulls `var2` up from below, `var_Δ5y` pulls it
down from above, and the canonical free estimate sits almost exactly where
these two pulls balance.

### 23b. The same profile for Germany: no interior minimum, a clean visual contrast

`scripts/ger_beta2_profile.py`, same full budget, `var2 ∈ {0.3, 0.5, 0.91,
1.5, 2.5}` (0.91 is GER's own canonical pin; GER's free/unrestricted
estimate sits at `var2 = 7.24`, far to the right of this grid). Result:

| var2 | OOS f | worst moment |
|---|---|---|
| 0.30 | 0.2912 | `var_log_earns` −53.2% |
| 0.50 | 0.1398 | `var_log_earns` −35.6% |
| 0.91 (KMV anchor, canonical) | 0.0354 | `var_log_earns` −14.7% |
| 1.50 | 0.0098 | `var_Δ5y` +6.2% |
| 2.50 | 0.0038 | `frac_d1_lt_20` −4.4% |

**Monotonically decreasing across the entire tested range, with no interior
minimum anywhere in `[0.3, 2.5]`** — a qualitatively different shape from
Poland's clean V with a well-defined floor at `0.509` (§23a). `var2=0.91`
reproduces the canonical GER pipeline's own OOS fit (`f ≈ 0.035`,
STEP2_LOG/§8b) almost exactly despite warm-starting from the free estimate
rather than KMV Table 3 — a useful independent cross-check that this profile
script's methodology matches the main pipeline's.

**Correction (2026-10-12): "no interior minimum" does not mean the curve is
flat or that the data are indifferent to `var2` — it is not, and they are
not.** The original wording here overstated the parallel with a genuinely
flat objective. Properly stated: **the data clearly reject low `var2`, and
clearly want it high — they just never say "high enough, stop."**
Quantified: pinning at KMV's own anchor (`var2=0.91`, the value actually
used for GER's canonical estimate) costs **14.0× in OOS fit** relative to
the free estimate evaluated on the identical panel (`f=0.0354` vs
`f=0.00252`, consistently computed — the free-estimate number used for this
ratio was itself corrected during this check, from an earlier
inconsistent-methodology evaluation of `0.0022`), and misses `var_log_earns`
by **−14.7%**. That is a real, substantial, single-moment-concentrated cost,
not a rounding error. The data's preference keeps strengthening without
leveling off all the way to the top of the tested grid (`var2=2.5`,
`f=0.0038`, still 1.5× worse than free estimation) — i.e. **the data require
`var2 ≳ 2.5`** to get anywhere near the best attainable fit, and the true
optimum (the free estimate, `var2=7.24`) sits further out still. What is
genuinely unidentified is only the **upper bound**: nothing in the data
rules out `var2` continuing to rise indefinitely as `β₂ → 0`, because the
objective keeps improving in that direction with no sign of turning around.
A lower bound is well pinned; a ceiling is not.

**Consequence: pinning `var2` at `0.9103` is a genuine modelling decision
that trades away real fit for a believable long-run dispersion (§8b's
calibration argument), not a value the SMM objective is indifferent to or
that the data themselves suggest.** This should be stated plainly wherever
the GER pin is defended — "the data don't care" and "we chose plausibility
over an unidentified free parameter" are different claims, and only the
second is accurate here.

**A hypothesis for why the data push `var2` ever higher, worth investigating
rather than asserting: a genuine trend in German earnings inequality over
2001–2010, not (or not only) a stochastic near-random-walk component.**
Checked directly: GER's raw `var(log annual earnings)` rises sharply and
almost monotonically from `0.559` (2001) to `0.666` (2010) — a 19% increase
over nine years — then flattens almost completely, `0.666` (2010) to
`0.674` (2016) (`output/pl_diagnostics/year_panel_ger.csv`). A stationary
jump-drift process has no way to represent a one-time, within-sample *rise*
in cross-sectional dispersion other than by inflating the persistent
component's ergodic variance — the stationary model's only lever for
"dispersion keeps growing" is a component that barely mean-reverts, i.e.
exactly the `β₂ → 0` pull this profile documents. If the true German data
generating process is non-stationary over this window (plausibly linked to
the well-documented rise in German wage inequality through the 2000s —
labour-market reforms, declining collective-bargaining coverage, and
low-wage-sector growth are the usual candidates in that literature, though
not verified here) rather than a stationary near-permanent shock, this would
be a textbook trend-mistaken-for-a-unit-root confound: the model absorbs a
deterministic or slow-moving trend by pushing a stochastic component toward
non-stationarity instead. **Not established here — a hypothesis for the
supervisor conversation and, if pursued, a candidate explanation for why
`β₂` resists identification specifically for Germany and not the other three
countries**, none of which show a comparably sharp early-sample level trend
(Poland's own `var_log` in fact *declines* over its sample, §18's diagnostics
addendum).

The combined PL-vs-GER picture (`output/pl_diagnostics/11_var2_profile.png`)
still makes the identification-strength contrast visible at a glance —
Poland's curve has a floor in the range that matters; Germany's keeps
falling — but the picture is "a real, one-sided pull with no ceiling," not
"indifference."

---

## 24. Poland's `var_log_earns` target measured the wrong variable — `std_perm` vs `std_loginc` (2026-10-19)

All §18–§23 Polish results (tag `notatka-promotor-v1`) used `targets/pol_grid.json`,
whose `var_log_earns = 0.5078` is the mean over years of `std_perm²`. FRA/GER/USA
use `std_log_inc²`. This section documents the problem, the fix and its
consequences. v1 outputs are left untouched; everything new goes to
`output/pol_v2/` and `targets/pol_grid_v2.json`.

### 24.1 The problem, confirmed

**What `std_perm` is (GRID Codebook, `grid-database.org/.../GRID_Codebook.pdf`,
"Residual log permanent income", p. 48).** The sample standard deviation of
*residual log permanent earnings*. Two things differ from annual log earnings,
not one:

1. **3-year average.** Permanent earnings = sum of annual earnings over
   `t−2 … t` divided by the number of those years above the minimum threshold,
   kept only if at least two years qualify. It is the average of *levels*,
   then logged. That is consistent with Guvenen–Pistaferri–Violante (NBER
   w30524, §2, "Variable definitions", item 3). The empty `std_perm` for
   2000–2001 in the file (no `t−2` year) is consistent with this.
2. **Age-residualized.** `r_i` is "the residual from the sex-specific
   age-dummy regression of log permanent earnings, estimated on sample [PE]
   in the relevant year". The 3-year average is regressed on age dummies
   within each year, so between-age differences in mean log earnings are
   removed.

It is also measured on a **different sample** ([PE]: ≥ 2 of 3 years above the
threshold, `nobs_perm` ≈ 4.1–4.5M) than annual log income ([M]: above the
threshold in year t, `nobs_loginc` ≈ 4.5–4.9M). [M] is the larger sample and
includes the weakly attached, low-earnings entrants that [PE] drops.

**What `std_loginc` is (Codebook p. 31).** The sample standard deviation of
the log of eligible annual earnings (above the threshold) on sample [M].
It is **raw**, not residualized: the residualized annual counterpart is a
separate statistic, `std_resid` (p. 40).

**The FRA/GER/USA `std_log_inc` is the same variable.** The July 2026 export
uses an older naming scheme, but it marks every residualized statistic with
`res_` (`std_res_1yr_log_chg`, …) and the level column without it
(`std_log_inc`), exactly as the October scheme does (`std_resid_*` vs
`std_loginc`). Sample sizes follow the same pattern: in FRA 1991,
`nobs_log_inc` = 302,370 > `nobs_res_1yr_log_chg` = 273,094. In PL 2010,
`nobs_loginc` = 4.74M > `nobs_resid_d1yr`. The level count is the larger [M]
cross-section, not the longitudinal sample. So `std_log_inc` (FRA/GER/USA) ≡
`std_loginc` (PL): **raw annual log earnings, not residualized.** This
verification is by naming convention and sample counts, since the old export
has no codebook of its own. It is not a byte-level check.

**What the model computes.** `moments_from_panel` takes
`var_log_earns = mean over panel columns of the cross-sectional var(log annual earnings)`.
That is the annual log level, **within each age** (one column = one age),
averaged over ages.

**Simulation (`scripts/pol_perm3_check.py` → `output/pol_v2/perm3_check.json`).**
v1 Polish parameters, 50k workers, seed 1234, the estimation panel. The
3-year average follows GRID's construction: average of levels over
`t−2..t`, then log. The threshold rule never binds in the model.

| model quantity (within-age, averaged over ages) | value |
|---|---|
| var(log annual earnings), all 36 columns (= `moments_from_panel`) | 0.5049 |
| var(log annual earnings), columns 3–36 | **0.5248** |
| var(log 3-year average), columns 3–36 | **0.4695** |
| ratio | **0.895** |

This reproduces the pre-check (0.470 vs 0.525, about 10% lower). **In the
model, 3-year averaging alone reduces the variance by about 10.5%.**

**In the data the gap is four times larger.** Same file
(`data/PL std_log_inc.csv`, whose `std_perm` column is identical to
`data/dane grid PL.csv` to the last digit, checked), mean of `std²` over years:

| window | `std_perm²` (v1 target) | `std_loginc²` | ratio |
|---|---|---|---|
| 2000–2020 (`std_perm` effectively 2002–2020) | 0.5078 | **0.7216** | 0.704 |
| 2002–2020 | 0.5078 | 0.7221 | 0.703 |
| 2001–2016 (the FRA/GER/USA window) | 0.5186 | 0.7386 | 0.702 |

In the data, `std_perm²` is 30% below `std_loginc²`; in the model, 3-year
averaging explains 10.5%. The other ~20 p.p. come from the two remaining
differences: age residualization (the between-age variance of mean log
earnings, which `std_perm` removes and `std_loginc` keeps) and the [PE] vs [M]
sample (weakly attached workers with low log earnings are in [M], not in
[PE]). This file does not separate the two, since it has no age breakdown and
no `std_resid`. **The correction to the target is therefore +42%
(0.508 → 0.722), not about +10%.** The v1 target was not "slightly too low";
it measured a much narrower concept.

**A separate issue, not caused by the Polish file.** Because
`moments_from_panel` computes the variance *within age* (per column), the
model moment conceptually corresponds to *age-residualized* log annual
earnings (GRID `std_resid`), not raw `std_loginc`. Raw `std_loginc` also
contains the between-age variance of mean log earnings, which the model, with
no deterministic age profile, does not generate. This applies identically to
**all four** countries (FRA/GER/USA have used raw `std_log_inc` from the
start, §10). So it does not break cross-country comparability. It does mean
`var_log_earns` is biased upward for all of them by the between-age variance.
It is a separate methodological choice from the one below (see 24.2).

### 24.2 Variant chosen: A

`data/PL std_log_inc.csv` (GRID Poland, Male, 25–55, 2000–2020, columns
`std_perm` and `std_loginc`) is in `data/`, so variant A applies.
`targets/pol_grid_v2.json` = `targets/pol_grid.json` with
`var_log_earns = mean_{2000..2020}(std_loginc²) = 0.72155`. The other seven
moments are copied unchanged byte for byte. The window 2000–2020 is the same
as the other seven Polish moments (§18). Results are insensitive to choosing
2002–2020 instead (0.7221).

### 24.3 Between-age variance: `std_resid²` vs `std_log_inc²`, all four countries

**Status of v2 at this point.** The v2 estimate (`output/pol_v2/`, pt 3) is
done: DE canonical, OOS f = 0.0073, var₂ = 1.13. The var₂ profile on v2
targets (pt 4) was **stopped** while its first point was still running, and
nothing from it was kept. If the level target changes again (below), the
profile and possibly the pt 3 estimate are stale. So pt 4 waits for this
decision.

**Data.** `data/std resid.csv`: GRID `std_resid` (Codebook p. 40: standard
deviation of the residual from a year-specific regression of log annual
earnings on age dummies, by sex). Male, 25–55, all four countries. The
file codes Germany as `DEU`. `std_log_inc` comes from the same sources as
the targets. Windows are the target windows: FRA/GER/USA 2001–2016, PL
2000–2020. Between-age variance = `mean(std_log_inc²) − mean(std_resid²)`.
Script: `scripts/resid_vs_loginc.py` → `output/resid_vs_loginc.json`.

| country | window | mean `std_resid²` | mean `std_log_inc²` | between-age | share of `std_log_inc²` | current `var_log_earns` target |
|---|---|---|---|---|---|---|
| POL | 2000–2020 | **0.6912** | 0.7216 | 0.0304 | 4.2% | 0.7216 (v2); v1 `std_perm²` 0.5078 |
| FRA | 2001–2016 | **0.4442** | 0.4762 | 0.0321 | 6.7% | 0.4762 |
| GER | 2001–2016 | **0.5890** | 0.6452 | 0.0561 | 8.7% | 0.6452 |
| USA | 2001–2016 | **0.8885** | 0.9557 | 0.0672 | 7.0% | 0.9557 |

- **The current targets are exactly the raw variance.** For FRA/GER/USA and
  PL v2, `var_log_earns` = mean `std_log_inc²` to the last digit. This also
  confirms that the July export's `std_log_inc` is the raw variable: it
  differs from `std_resid` in every year and every country, and is always
  larger.
- **The between-age component is small and stable.** By year it ranges
  0.020–0.042 (PL), 0.028–0.037 (FRA), 0.040–0.061 (GER) and 0.048–0.081
  (USA). In the USA it rises over the window (0.048 in 2001, 0.080 in 2013).
  The country ranking does not change: PL < FRA < GER < USA in absolute terms.
- **It does not explain the PL v1 → v2 jump.** Between-age variance in Poland
  is only 0.030. Most of the 0.21 gap between `std_perm²` and `std_loginc²`
  (24.1) therefore comes from 3-year averaging (the model says about 0.06)
  and from the [PE] vs [M] sample (the remaining ~0.12), not from
  residualization.
- **Same sample.** `nobs_resid` = `nobs_log_inc` exactly for PL/FRA/USA in
  every year, so `std_resid` and `std_log_inc` are on the same sample [M].
  The difference between them is pure age residualization.
  **Exception: GER.** The July export's `nobs_log_inc` is exactly 3.997× the
  new export's `nobs_resid`, at a constant ratio (0.2502) in every year. A
  constant ratio points to a different count *reporting* convention in one of
  the two exports (e.g. weighting / subsample of the IAB data), not to a
  different selection. That would change the ratio from year to year. Not
  confirmed: the clean check is to download `std_loginc` for GER from the same
  new export.

**What this means for the comparison with the model.** The model moment
(`moments_from_panel`) is the within-age variance averaged over ages, so
conceptually it matches `std_resid²`, not `std_log_inc²`. Raw targets
overstate what the model should match by 4–9%, depending on the country.

**A second mismatch on the model side, in the opposite direction, which
must be decided together.** The simulated panel covers **36 years**
(≈ ages 25–60, §2, chosen to replicate KMV/GKOS). GRID covers ages
**25–55**, i.e. 31 years. Within-age variance grows with age, so the
extra 5 oldest years raise the model's `var_log_earns`. Simulated at each
country's canonical parameters (50k workers, seed 1234), mean variance over
the first 31 columns vs all 36:

| parameters | all 36 years | first 31 years (25–55) | difference |
|---|---|---|---|
| PL v2 | 0.720 | 0.670 | −6.9% |
| FRA | 0.476 | 0.441 | −7.3% |
| GER | 0.551 | 0.510 | −7.3% |
| USA | 0.803 | 0.770 | −4.1% |

The two biases have similar magnitude and **opposite signs**. Raw data add
between-age variance of +4–9%. The model adds ages 56–60, worth +4–7%. The
current setup (raw target, 36-year panel) is therefore roughly
self-cancelling, by coincidence rather than by design. **Switching to
`std_resid` alone, without trimming the panel to 31 years (or computing
moments only on the 25–55 columns), would replace one bias with another of
similar size.** Consistent options:
1. **`std_resid` + model moments on ages 25–55.** Fully consistent on both
   sides. Requires re-estimating all four countries, because the panel
   change also affects `var_d1`/`var_d5`/kurtosis (by a smaller amount, not
   measured here) and the targets change for all of them.
2. **Status quo (raw `std_log_inc` + 36 years) as a documented approximation,**
   with the table above as justification that the net bias is a few percent.
   PL v2 stays as is.

Recommendation and decision: in 24.4, after discussion.

### 24.4 Decision: option 2, raw `std_log_inc` + 36-year panel for all countries (2026-10-19)

**Decision.** Keep the level-dispersion target as raw `std_log_inc²` (GRID
annual log earnings, sample [M]) and keep the 36-year simulation panel
(§2) for **all four countries**. Poland uses `targets/pol_grid_v2.json`
(24.2). FRA/GER/USA targets are unchanged.

**Rationale.**
1. **The two biases from 24.3 have similar magnitude and opposite signs, net
   a few percent.** Raw data contain between-age variance that the model
   lacks: +4.2% (PL), +6.7% (FRA), +8.7% (GER), +7.0% (USA) of the target.
   The 36-year panel adds ages 56–60, which GRID does not have: +6.9% / +7.3%
   / +7.3% / +4.1% of the model moment. The net mismatch is about −3 to +3%
   depending on the country. That is within the Monte Carlo noise of the SMM
   objective at this scale (§9) and smaller than the fit error on the worst
   moment in every country.
2. **Cross-country consistency is preserved.** All four countries use the
   same GRID variable for the target and the same model moment. That was the
   validity criterion from the start (§11/§12). Fixing just one of the two
   biases, i.e. `std_resid` without trimming the panel, would *break* this
   balance.
3. **Option 1 means re-estimating all four countries,** and it changes the
   change-moment definitions too (the 25–55 panel affects `var_d1`/`var_d5`
   and kurtosis). That is a large cost for a gain of a few percent on one
   moment. It belongs as a robustness test, not as the baseline.

**Option 1** (`std_resid` + model moments on ages 25–55) goes to "Open
decisions" as a robustness test to discuss with the supervisor.

### 24.5 Report: Poland v1 vs v2 (2026-10-20)

v1 = `targets/pol_grid.json` (`var_log_earns` = `std_perm²`, 0.5078), outputs
`output/pol/`, `output/pol_var2_*` (tag `notatka-promotor-v1`).
v2 = `targets/pol_grid_v2.json` (`var_log_earns` = `std_loginc²`, 0.7216),
outputs `output/pol_v2/`. Same pipeline (`close_step1.py`: DE + TikTak,
50k workers, DE maxiter 100 / popsize 12, NM 400, warm start at KMV Table 3;
OOS on 100k workers, seed 987654). Profile `pol_beta2_profile.py --version v2`,
same budget, warm start at the v2 canonical.

**Parameters (canonical = DE in both versions)**

| | v1 | v2 | change |
|---|---|---|---|
| λ₁ | 0.1623 | 0.1623 | 0% |
| β₁ | 2.846 | 2.832 | −0.5% |
| σ₁ | 1.979 | 2.001 | +1.1% |
| λ₂ | 0.00877 | 0.00987 | +12.5% |
| β₂ | 0.01690 | **0.00638** | **−62%** |
| σ₂ | 1.401 | 1.209 | −13.7% |
| ergodic sd z₂ (persistent) | 0.714 | **1.064** | +49% |
| ergodic sd z₁+z₂ | 0.788 | 1.117 | +42% |
| half-life, transitory | 0.061 yr (≈ 3 weeks) | 0.061 yr | — |
| half-life, persistent | 10.3 yr | **27.2 yr** | ×2.65 |
| OOS f (DE / TikTak) | 0.00367 / 0.00370 | **0.00728** / 0.00818 | ×2 |
| DE vs TikTak, persistent component (max rel. diff) | 1.3% | 12.3% | still AGREE (< 25%) |
| parameter at a bound | none | none | |

The transitory component is unchanged to within 1%. The entire +0.21 of the
level target was absorbed by the persistent component, almost entirely
through slower mean reversion (β₂ fell by 62%), with partial compensation
from a smaller σ₂ and slightly more frequent shocks. This is exactly the
direction the v1 profile predicted (§23: `var_log_earns` disciplines var₂
from below).

**Fit to the 8 moments (OOS, rel. error model vs data)**

| moment | data v1 | model v1 | err v1 | data v2 | model v2 | err v2 |
|---|---|---|---|---|---|---|
| var log earnings | 0.5078 | 0.5057 | −0.4% | **0.7216** | 0.7211 | −0.1% |
| var Δ1y | 0.2595 | 0.2595 | 0.0% | 0.2595 | 0.2627 | +1.2% |
| var Δ5y | 0.4765 | 0.4835 | +1.5% | 0.4765 | 0.4771 | +0.1% |
| kurt Δ1y | 15.24 | 14.80 | −2.9% | 15.24 | 14.52 | **−4.7%** |
| kurt Δ5y | 9.878 | 9.979 | +1.0% | 9.878 | 9.982 | +1.1% |
| frac \|Δ1y\| < 10% | 0.4499 | 0.4710 | +4.7% | 0.4499 | 0.4798 | **+6.6%** |
| frac \|Δ1y\| < 20% | 0.6694 | 0.6577 | −1.7% | 0.6694 | 0.6575 | −1.8% |
| frac \|Δ1y\| < 50% | 0.8431 | 0.8395 | −0.4% | 0.8431 | 0.8363 | −0.8% |

The doubling of OOS f comes almost entirely from two moments: `frac_d1_lt_10`
(contribution 0.0022 → 0.0044) and `kurt_d1` (0.0008 → 0.0022). This misfit
is **constant along the whole v2 profile** (frac +6.5–6.6%, kurt −4.0 to −4.7%
for var₂ from 1.0 to 1.6, below), so it does not come from var₂. It comes from
the new target as a whole: a persistent component with more frequent but
smaller shocks (λ₂ ↑, σ₂ ↓) adds small annual changes, so the model
overstates the share of |Δ1y| < 10% and lowers the kurtosis. With the
transitory parameters essentially fixed, the optimiser cannot compensate.
In absolute terms the fit is still good: no moment is off by more than 6.6%.

**var₂ profile (`output/pol_v2/var2_profile_summary.json`, plot
`output/pol_v2/11_var2_profile_v2.png`, with the GER curve and v1 faded)**

| var₂ | OOS f | vs min | worst moment (excluding the constant frac10/kurt1 misfit) |
|---|---|---|---|
| 0.6 | 0.0540 | ×7.4 | `var_log_earns` −21.2% |
| 0.8 | 0.0167 | ×2.3 | `var_log_earns` −8.8% (and `var_d5` +5.9%) |
| 1.0 | 0.0080 | +9% | — (`var_log_earns` −2.1%, `var_d5` +2.2%) |
| **1.132 (control)** | **0.00728** | 0 | identical to the free estimate (f to 5 digits, params ≤ 0.1%) |
| 1.3 | 0.0079 | +9% | — (`var_d5` −2.0%) |
| 1.6 | 0.0099 | +36% | `var_d5` −4.4% (and `var_d1` +3.1%) |

Noise band (+20% of the minimum, §9): 0.0087. Points 1.0, 1.132 and 1.3 are
inside it. 0.8 and 1.6 are outside. The control point reproduces the free
estimate, so the pin does not bind at the optimum and the profile
methodology holds (as in §23 for v1).

**Worst moment on both sides of the profile.** Left (var₂ too low):
`var_log_earns`, the same as in v1 (v1: −26% at 0.3, −13% at 0.4). Right
(var₂ too high): `var_d5` (−2.0% at 1.3, −4.4% at 1.6), the same as in v1
(v1: −8.9% at 0.7, −13.5% at 0.91). So the identifying structure has not
changed: the level pulls var₂ up, the 5-year change pulls it down. Only the
point where the two balance has shifted. In v2 the right slope is much
**flatter** than in v1: going ×1.4 above the optimum costs +36% of f in v2,
and in v1 (0.509 → 0.7) it cost ×5.5.

**Answers to the three questions**

1. **Does Poland still have the lowest long-run dispersion of the four
   countries? No.** Ergodic sd of the persistent component: PL v2 **1.064**
   vs FRA 0.903, USA 0.954 (pinned; free 1.22), GER 0.954 (pinned; free
   2.69). Combined: PL v2 1.117 vs FRA 0.947, GER 0.979, USA 1.024. Among the
   canonical estimates Poland now has the **highest** dispersion. It is still
   below the free USA and GER estimates, but those are not canonical (§8b, S5).
   In the data the ordering was always clear and v1 hid it: annual
   `var_log_earns` PL 0.722 > GER 0.645 > FRA 0.476, behind only the USA
   (0.956). The v1 conclusion that Poland is the most mean-reverting and has
   the lowest dispersion (§18) was an artefact of the `std_perm` target and
   should be withdrawn.
2. **Do the data still clearly prefer a value below the KMV anchor (0.91)?
   No. The direction is reversed.** The free optimum is var₂ = 1.13, above
   0.91. The anchor itself was not computed on the v2 grid. Log-linear
   interpolation between 0.8 and 1.0 gives f(0.91) ≈ 0.011 (+50% over the
   minimum, outside the noise band), so the data mildly reject 0.91 *from
   below* (`var_log_earns` −5% would be the worst moment). In v1, 0.91 was
   rejected strongly from above (f ×8.4, `var_d5` −13.5%). If an exact
   number is needed for the write-up, one run at `--var2 0.9103` (~35 min)
   settles it.
3. **Does the profile still have a clear floor? Yes, but wider and
   asymmetric.** Two-sided: 0.8 (×2.3) and 1.6 (+36%) are outside the noise
   band, and the floor is roughly var₂ ∈ [≈0.95, ≈1.4], i.e. ergodic sd
   ≈ 0.97–1.18. Compared with v1, the floor is relatively wider (in v1 only
   the canonical point was in the band; 0.4 and 0.7 were already ×5.5), and
   the right slope is shallow. This is still the opposite of Germany, whose
   curve has no minimum and keeps falling up to var₂ = 7.2. The PL/GER
   contrast on the plot holds, though it is weaker than in v1.

**Does this change the step-2 numbers?** Yes. See 24.6 (`swap_death`). v2
`tables/10_economic_params.csv` has a new `POL_V2` column (the `POL` column =
v1, left for comparison).

### 24.6 Two-asset trigger on the v2 chain: Poland does NOT calibrate (2026-10-20)

`hank_step2/swap_death.py pol_v2`: the same procedure as §20 (one-asset ARS +
β heterogeneity + KMV death, joint calibration of (β_hi, dβ, ω) to A = 20,
impact labor MPC = 0.20, SCF Lorenz), on the chain `output/pol_v2/`.
Result: `hank_step2/results/swap_death_pol_v2.json` (+ `.log`; a `pol_v2` key
merged into `swap_death_summary.json`, other keys byte-identical).

| | target | **POL v2** | POL v1 (§20) | FRA | GER | USA |
|---|---|---|---|---|---|---|
| converged | — | **no**, cost 0.0087 | yes, resid ~1e-12 | yes | yes | no |
| MPC (labor, impact) | 0.20 | **0.071** | 0.200 | 0.200 | 0.200 | 0.024–0.029 |
| frac at borrowing constraint | — | **0.063** | 0.226 | 0.225 | 0.250 | ~0.01 |
| dβ (heterogeneity spread) | — | 0.150 | 0.125 | 0.093 | 0.067 | →0 |
| β_lo | — | 0.843 | 0.873 | 0.905 | 0.932 | — |
| ω (share of patient type) | — | **0.934** | 0.691 | 0.699 | 0.666 | 0.592 |
| decomp yr1 direct/indirect | — | **35/65** | 22.6/77.4 | 25.5/74.5 | 24/76 | degenerate |
| monetary dY impact (base 2.41) | — | **2.06 (−14.5%)** | 2.30 (−4.5%) | 2.29 | 2.33 | 1.81–1.85 |
| ergodic var(log e), chain | — | **1.381** | 0.721 | 1.040 | 1.044 | 1.191 |

**This is not a local minimum.** The solver stopped inside the box: dβ = 0.15
vs a limit of 0.30, ω = 0.93 vs 0.98. So the calibration step alone was
rerun from 6 more starts (`hank_step2/swap_death_multistart.py`,
`results/swap_death_pol_v2_multistart.json`). They included the converged
solutions for PL v1 and for var₂ = 0.7 and the corners of the box. **All 9
starts (3 + 6) end on the same ridge:** β_hi = 0.9930, ω = 0.927–0.934,
dβ = 0.105–0.150 (flat in dβ), with an MPC residual of −0.129 to −0.142,
A = 20 matched to 0.003 and Lorenz to −0.02 to −0.03. In this model,
the v2 chain can deliver at most MPC ≈ 0.07 while keeping A and Lorenz.

**Mechanism (same as GRID-USA, STEP2_LOG §S4/§S5).** Ergodic variance of
income in the chain is 1.38. That is about double v1 and a third above
FRA/GER. The persistent component's half-life is 27 years. With such large
and durable income risk, precautionary saving pushes almost everyone away
from the borrowing constraint (6% vs 23% in v1). The calibration tries to
recreate a hand-to-mouth group with an impatient type, but to hit A = 20 and
Lorenz it must keep that type small (1 − ω = 7%). MPC does not get above
0.07. The monetary transmission is weaker (−14.5% vs the baseline) and the
indirect share drops to 65%. This is the USA failure mode, though milder
(MPC 0.07 vs 0.03).

**Note: the dependence on var₂ is non-monotonic.** On the v1 targets, the
var₂ = 0.91 chain (ergodic var 1.08) failed in the *opposite* direction:
dβ at the 0.30 limit, MPC 0.33, HtM 0.34 (§23). Polish chains therefore
calibrate in a window of ergodic variance roughly between 0.72 (v1, OK) and
0.95 (v1 var₂ = 0.7, OK) — v1 var₂ = 0.91 already overshoots, while v2 at
1.38 undershoots. These are chains with different targets, so they are not
one comparative static. The v2 profile floor, however, covers var₂ ∈
[≈0.95, ≈1.4]. **Untested: whether a chain from the lower part of the v2
floor (`output/pol_v2/var2_1.0/`, OOS f +9% over the minimum, inside the
noise band) calibrates.** One `swap_death.py` run (~2 h) answers this, and it
matters for the two-asset decision (below).

**Consequence.** The §20 conclusion ("Poland calibrates cleanly, one-asset +
death is enough") **applies to the v1 targets and does not hold for v2.**
The two-asset question is back in "Open decisions".

### 24.7 Do the bootstrap, recovery, nested specs and sub-period split need rerunning for Poland v2? (assessment, not run)

**Short answer: yes, all four.** Marking them as "for the previous version"
is enough for the supervisor note *for now*, but none of them can carry over
to v2 as a result. The parameter shift is not small, and in two cases v1
itself depended on `std_perm`.

| exercise | v1 result (§) | why it does not carry over to v2 | priority / cost |
|---|---|---|---|
| **Sub-period split pre/post 2009** (§20 Q2, `pol_regime_split.py`) | the two-asset answer is the same in both sub-windows | `targets/pol_grid_{pre,post}2009.json` are built from `dane grid PL.csv`, so `var_log_earns` = `std_perm²` (0.541 / 0.488). **The same error as in the full window.** The conclusion ("calibrates in both") also cannot hold, since the full sample now does not calibrate (24.6). | **1**, ~2 × (1 h estimation + 2 h swap_death) |
| **Nested specs** (§22, `nested_specs.py`) | V3 (β₂ = 0, random walk) is clearly worse: OOS f 0.0130 vs 0.0037 (×3.5) | v2 moves β₂ toward zero (half-life 10 → 27 years). The key conclusion that the persistent component must mean-revert may weaken or flip. **This is the test most sensitive to exactly this change.** | **2**, ~4 estimations ≈ 2–3 h |
| **Recovery** (§21, `recovery_exercise.py`) | β₂ recovered with mean bias −25% and a wide range (min −99%) — identification already weak in v1 | recovery tests identification **at a given point**. The v2 point is in a different region (β₂ ×0.38, λ₂ +12%, σ₂ −14%) and moves toward the weakly identified German region (§8b). A v1 result says nothing about v2, and β₂ identification is likely worse there. | **3**, as in §21 |
| **Bootstrap SEs** (§19, `bootstrap_se.py`) | β₂ 95% CI [0.0150, 0.0193], λ₂ [0.0082, 0.0096], σ₂ [1.383, 1.452] | **The v2 point lies far outside the v1 CIs** for all three persistent-component parameters (β₂ = 0.0064, λ₂ = 0.0099, σ₂ = 1.209), so the v1 SEs describe a different estimate. They also need a new year panel: `output/pl_diagnostics/year_panel_pol.csv` has `var_log` from `std_perm` (empty for 2000–2001). The panel must be rebuilt from `std_loginc` before the bootstrap. | **4**, the most expensive (200 reps on a reduced budget, as in §19) |

The transitory component (λ₁, β₁, σ₁) has not moved (≤ 1%), so the v1
conclusions about it (bootstrap SEs, recovery of β₁/σ₁) are **probably** still
valid. But they cannot be reported without a rerun, because the parameters
are estimated jointly.

**Recommended order:** (0) `swap_death` on `output/pol_v2/var2_1.0/` (24.6,
~2 h), because it decides whether Poland goes to two-asset, and that frames
everything else. Then (1) sub-period split, (2) nested specs, (3) recovery,
(4) bootstrap. In the supervisor note: the §19–§23 Polish results are marked
as "v1 targets (`std_perm`), superseded — v2 rerun in progress", with §24 as
the explanation.

### 24.8 Sub-period split pre/post 2009 on v2 targets (2026-10-20)

**Targets.** `targets/pol_grid_{pre,post}2009_v2.json` = the v1 sub-period
files with `var_log_earns` replaced by `mean(std_loginc²)` over the sub-window
(`data/PL std_log_inc.csv`): pre 2000–2008 **0.7735** (v1 `std_perm²`:
0.5411, effectively 2002–2008), post 2009–2020 **0.6826** (v1: 0.4884). The
other 7 moments are copied unchanged. Estimation is exactly as in §20
(`pol_regime_split.py --version v2`: DE, 50k workers, maxiter 100 / popsize
12, NM 400, start at KMV Table 3), with **no** `swap_death`. Added: OOS f on
100k workers, seed 987654. Outputs: `output/pol_{pre,post}2009_v2/`,
`output/pol_regime_split_v2.json`.

| | full v1 | pre v1 | post v1 | **full v2** | **pre v2** | **post v2** |
|---|---|---|---|---|---|---|
| `var_log_earns` target | 0.508 | 0.541 | 0.488 | 0.722 | 0.773 | 0.683 |
| λ₁ | 0.162 | 0.177 | 0.153 | 0.162 | 0.177 | 0.153 |
| β₁ | 2.85 | 2.86 | 2.82 | 2.83 | 2.89 | 2.83 |
| σ₁ | 1.98 | 1.99 | 1.98 | 2.00 | 2.02 | 1.99 |
| λ₂ | 0.00877 | 0.00960 | 0.00775 | 0.00987 | 0.01138 | 0.00853 |
| β₂ | 0.0169 | 0.0178 | 0.0130 | 0.0064 | 0.0061 | **0.0044** |
| σ₂ | 1.40 | 1.40 | 1.33 | 1.21 | 1.15 | 1.18 |
| ergodic sd z₂ | 0.714 | 0.726 | 0.728 | 1.064 | 1.113 | 1.166 |
| half-life z₂ (yrs) | 10.3 | 9.8 | 13.4 | 27.2 | 28.6 | **39.6** |
| in-sample f | 0.0041 | 0.0046 | 0.0032 | 0.0078 | 0.0085 | 0.0048 |
| OOS f | 0.0037 | — | — | 0.0073 | 0.0082 | 0.0045 |

Fit (OOS) in the sub-periods: worst moments are `frac_d1_lt_10` (+6.7% pre /
+5.5% post) and `kurt_d1` (−5.6% / −3.3%), the same pattern as the full v2
window (§24.5). `var_log_earns` is matched to ≤ 0.5%. No parameter is at a
bound.

**Readings.**
1. **The transitory component is the same in all six estimates** (β₁
   2.82–2.89, σ₁ 1.98–2.02). The only thing that changes between sub-periods
   is λ₁ (0.177 pre vs 0.153 post), and identically in v1 and v2. So the
   2009 regime break (kurtosis / small-change share, §20) shows up as less
   frequent transitory shocks, regardless of the level target.
2. **The level change v1 → v2 does the same thing in each sub-period as in
   the full window**: β₂ falls about 3× and the persistent ergodic sd rises
   from about 0.72 to 1.1–1.17. The v2 conclusion (24.5: Poland's dispersion
   is high, not low) holds in both sub-periods.
3. **Post-2009 v2 has the slowest mean reversion (half-life 40 years) despite
   a *lower* level target (0.683 vs 0.773 pre).** This is the same mechanism
   as in the profile (24.5): the level pulls var₂ up and `var_d5` pulls it
   down. Post-2009 `var_d5` is much lower (0.424 vs 0.517 pre), so the
   discipline from above is weaker and β₂ falls further. The 12-year window
   is also short for identifying a ~40-year half-life (cf. GER, §8b). The
   post-2009 β₂ figure should be treated as weakly identified. It is not
   evidence of a trend.
4. **Comparison with §20 Q2.** In v1, both sub-periods gave practically the
   same persistent sd (0.726 / 0.728) as the full window (0.714). In v2 both
   sub-periods give it **higher** than the full window (1.11 / 1.17 vs 1.06),
   and the spread between them is larger (half-life 29 vs 40 years). So the
   persistent component is less stable across sub-periods than v1 suggested.
   **Not checked here (as instructed): whether either sub-period calibrates
   in the one-asset model.** Given ergodic sd ≥ 1.11 (≥ the full v2 window,
   which does not calibrate, 24.6), it is very unlikely. The §20 Q2 answer
   ("calibrates in both sub-windows") does not carry over to v2.

### 24.9 Nested specs for Poland v2 (2026-10-20)

`nested_specs.py --version v2 --countries POL`. Same design and budget as
§22 (DE, 50k workers, maxiter 100 / popsize 12, NM 400, start at KMV
Table 3; OOS 100k workers, seed 987654). Reference: the free v2 estimate
(`output/pol_v2/params.json`). Outputs `output/nested_specs_v2/` (v1 in
`output/nested_specs/` untouched).

| variant | free params | OOS f v1 | ×full v1 | **OOS f v2** | **×full v2** | worst moment v2 |
|---|---|---|---|---|---|---|
| full (6) | 6 | 0.00367 | 1 | **0.00728** | 1 | `frac_d1_lt_10` +6.6% |
| V1 single component | 3 | 0.308 | ×84 | 0.481 | ×66 | `var_log_earns` −61% |
| V2 deterministic arrival (AR(1)) | 2 | 2.19 | ×596 | 2.36 | ×324 | `kurt_d1` −80% |
| **V3 β₂ = 0 (random walk)** | 5 | 0.0130 | **×3.5** | **0.0224** | **×3.1** | **`var_d5` −10.9%** |

**Key question: is V3 (β₂ = 0) still clearly worse than the full model?
Yes.** OOS f 0.0224 vs 0.0073 is ×3.1, far outside the noise band (+20% of
the minimum, 0.0087). Relative to v1 (×3.5) the penalty is barely smaller,
even though v2 moved β₂ toward zero (half-life 27 years). Mean reversion of
the persistent component therefore stays necessary.

**What changed in V3's structure, and why it makes the v2 test stronger.**
In v1 the V3 optimum was a degenerate corner: λ₂ at the 0.1 bound, σ₂ = 0.21,
and the "transitory" component took over persistence (β₁ = 0.50, λ₁ = 0.058).
It was a different model, not a random walk added to the same transitory
component. In v2, V3 has a sensible structure. The transitory component is
almost the same as the full one (λ₁ 0.153, β₁ 2.60, σ₁ 2.08). The
permanent component has λ₂ = 0.0133 (a shock every ~19 years) and σ₂ = 0.81.
So v2 compares like with like: the same transitory structure, with and
without β₂ > 0. **It is rejected by `var_d5`** (−10.9%, plus `var_d1` +6.2%,
`var_log_earns` +4.3%). A random walk cannot reproduce the average level
(0.72, an average over 36 years in which variance grows linearly with age)
and the 5-year change variance at the same time. Fitting the level forces a
small σ₂, and then the permanent part contributes too little to Δ5y. This is
the same discipline from above as on the right side of the profile (24.5).

**V1/V2:** still catastrophic. They are slightly worse in absolute terms
than in v1 (the higher level target is harder to hit with a single component:
`var_log_earns` −61% vs −44% in v1). Their ratio to the full model is lower
only because the denominator doubled. No change in conclusions: both
components and Poisson arrival are necessary.

### 24.10 Recovery for Poland v2 (2026-10-20)

`recovery_exercise.py --version v2 --countries POL`. The same scheme as §21:
the true process is the free v2 estimate (`output/pol_v2/params.json`),
2 generating seeds (101, 202) × 2 synthetic panel sizes (50k, 300k).
Re-estimation uses DE at the full budget, starting from the **neutral** KMV
Table 3 values (not from the truth). Outputs `output/recovery_v2/` (v1 in
`output/recovery/` untouched).

| rep (seed, panel) | β₂ reldiff v1 | β₂ reldiff **v2** | ergodic sd z₂ v1 (true 0.714) | ergodic sd z₂ **v2** (true 1.064) |
|---|---|---|---|---|
| 101, 50k | −0.9% | −1.3% | 0.713 | 1.065 |
| 101, 300k | −0.4% | −2.0% | 0.714 | 1.069 |
| 202, 50k | −0.2% | −1.6% | 0.713 | 1.068 |
| 202, 300k | **−99.4%** (failure: λ₂ ×11, a different basin) | −2.3% | **4.70** | 1.074 |

Mean over all parameters, v2: λ₁ −0.1%, β₁ −0.5%, σ₁ 0.0%, λ₂ −1.1%,
**β₂ −1.8% (sd 0.4%)**, σ₂ +0.1%; ergodic sd z₂ +0.5%. β₂ by panel size:
−1.5% (50k) vs −2.1% (300k), so it does not improve with sample size.

**Readings.**
1. **At the v2 point, β₂ is well identified, in the USA/FRA class (§21: −1.8%
   / −4.0%), not the GER class (−82%).** All 4 reps fall in one basin, with a
   small systematic bias of about −2% that does not shrink with the panel.
   That is an ordinary small-sample / simulation bias, not an
   identification problem. The fear in 24.7 ("v2 moves toward the German,
   weakly identified region") **is not confirmed**: a 27-year half-life is
   still identifiable from 8 moments on a 36-year panel. Germany's (§8b) is
   ~250 years.
2. **v1's bimodality (1/4 catastrophic reps) does not recur.** Four draws
   do not prove that failure is impossible: 0/4 is consistent with a
   failure rate of up to ~50% at 95% confidence. But it is consistent with
   v2's minimum being easier for DE to find from a neutral start than v1's.
3. **How this fits with the profile's shallow right slope (24.5).** Recovery
   asks whether the moments map back to the parameters when the *only* noise
   is simulation noise. The profile asks how much the fit to the *real* data
   worsens as var₂ moves. Both can be true at once: the parameters are
   point-identified (recovery), but the real data have a wide floor for var₂
   (1.0–1.3 in the noise band), i.e. moderate *statistical* precision. The
   bootstrap (24.11, if completed) gives the sampling-uncertainty part.

### 24.11 Bootstrap SEs for Poland v2 (2026-10-20)

**Panel.** `output/pl_diagnostics/year_panel_pol_v2.csv` is the v1 panel with
`var_log` = `std_loginc²` (and `nobs_level` = `nobs_loginc`) for each year
2000–2020. The v1 file is untouched. Check: the unresampled mean reproduces
`targets/pol_grid_v2.json` **to the bit** for `var_log_earns`, `var_d1/d5`
and `kurt_d1/d5`. The three small-change fractions in the panel (mean over
years of per-year fractions) differ from the target (fractions from averaged
percentiles) in the third digit (e.g. `frac_d1_lt_10` 0.4558 vs 0.4499).
**This is identical in v1** (checked: the same discrepancy in
`year_panel_pol.csv`), so it is a property of the §19 construction, not a new
error. It means the bootstrap reps are centred on fractions about 1% higher
than the point-estimate target. Not corrected, to keep the design "as in
§19".

**Run.** `bootstrap_se.py --version v2 --countries POL`, same budget as §19
(reduced: 5k workers, DE 10×5, NM 50, warm start at `output/pol_v2/params.json`).
Budget validation (3 resamples, reduced vs reference): OK. The largest mean
deviation is β₂ at −3.7% (sd 2.5%), below the 15% threshold. 200 reps,
1418 s. Outputs `output/bootstrap_v2/` (v1 in `output/bootstrap/` untouched).
The version is passed through `KMV_BOOTSTRAP_VERSION`, because pool workers
on Windows re-import the module and changes made in `main()` would not reach
them.

| | λ₁ | β₁ | σ₁ | λ₂ | β₂ | σ₂ |
|---|---|---|---|---|---|---|
| v2 point estimate | 0.1623 | 2.832 | 2.001 | 0.00987 | 0.00638 | 1.209 |
| v2 bootstrap mean | 0.1591 | 2.822 | 2.007 | 0.00981 | 0.00622 | 1.217 |
| v2 SE | 0.0033 | 0.077 | 0.010 | 0.00030 | 0.00060 | 0.012 |
| v2 95% CI | [0.153, 0.165] | [2.69, 2.91] | [1.989, 2.026] | [0.00938, 0.01045] | [0.00546, 0.00805] | [1.203, 1.251] |
| v2 rel. SE | 2.0% | 2.7% | 0.5% | 3.0% | **9.4%** | 1.0% |
| v1 rel. SE (§19) | 2.0% | 1.9% | 0.6% | 4.0% | 7.0% | 1.1% |

Derived quantities (per rep):

| | v1 | **v2** |
|---|---|---|
| var₂ 95% CI | [0.501, 0.540] | **[1.005, 1.276]** (mean 1.174) |
| ergodic sd z₂ 95% CI | [0.708, 0.735] | **[1.002, 1.130]** |
| half-life z₂ 95% CI (yrs) | [9.0, 11.6] | [21.5, 31.8] |

**Readings.**
1. β₂ is the least precise parameter, as in v1, and a little less precise
   than in v1 (relative SE 9.4% vs 7.0%; the CI is asymmetric, with the upper
   tail toward faster mean reversion). The other parameters have relative SEs
   of 0.5–3%.
2. **The var₂ CI [1.00, 1.28] excludes both the v1 value (0.51) and the KMV
   anchor (0.91).** That is firmer evidence than the profile interpolation
   (24.5, Q2) that the data reject 0.91. **Interpretation caveats, as in
   §19:** (i) this is only year-to-year variability of the GRID statistics,
   not full sampling uncertainty; (ii) every rep is warm-started at the point
   estimate, so the CI is conditional on that basin. §24.10 suggests the
   basin risk is small in v2 (4/4 recovery reps in one basin).
3. The bootstrap CI for var₂ is narrower than the profile's noise band
   (1.0–1.3, 24.5) and lies inside it, so the two measures agree. The
   `var2_1.0` chain (24.6 / §25, a candidate for the one-asset route) sits
   at the **lower edge** of the bootstrap CI.
4. The sub-period point estimates (24.8: var₂ 1.24 pre / 1.36 post) fall in
   the upper tail of this CI or above it. The difference between sub-periods
   is of the same order as the year-to-year variability. Again, this points
   to weak precision on the right side of the profile, not to a structural
   break.

---

## 25. Newborns' income distribution in HANK: ergodic draw vs start at z = 0 (2026-10-20)

**Question.** In `hank_step2/death.py`, newborns keep their (β, e) state, so
in the cross-section they draw `e` from the chain's ergodic distribution. For
PL v2 that distribution has var(log e) = 1.38, while the target is 0.72
(§24.6). The estimation panel (`simulate.py`) instead starts every worker at
z = 0. Would starting newborns at z = 0 close this gap? Is it consistent with
KMV? **Nothing in `death.py` has been changed.**

**What KMV do.** KMV (2018, AER), §II.A: households die at Poisson rate ζ
"and upon death give birth to an offspring with zero wealth and labor
productivity equal to a random draw from its ergodic distribution". Footnote
15 adds that death is there to generate enough households with zero
*illiquid* wealth, and that the stationary distribution exists even with
ζ = 0. ζ = 1/180 per quarter, average lifetime 45 years (Table 6). So
**`death.py` reproduces KMV exactly.** Their code was not checked
separately; the text is unambiguous. KMV themselves live with the same gap:
the process is estimated on a lifecycle panel from z = 0 (our `simulate.py`
reproduces their Table 2 that way), and the HANK cross-section is ergodic.
For KMV's own chain, ergodic var = 1.07 (STEP2_LOG S5; 1.15 after our
re-discretization of Table 3) vs the target of 0.70.

**Computation** (`scripts/newborn_entry_variance.py` →
`output/newborn_entry_variance.json`). Each country's exported chain (rebuilt
with the same `discretize_process` call, checked against
`income_process_zgrid.txt`). The components are independent, so the
cross-section is D_j = ζ·e₀(I − (1−ζ)P_j)⁻¹ with e₀ = point mass at z_j = 0.
The figures are var(log e) of **quarterly** productivity in the cross-section:

| chain | ergodic (now) | perpetual youth, start z = 0 | same, ages truncated at 36 yrs | uniform ages 0–36 yrs, z = 0 | target `var_log_earns` (annual) | ergodic / target |
|---|---|---|---|---|---|---|
| **POL v2** | **1.381** | **0.982** | 0.713 | 0.781 | 0.722 | 1.91 |
| POL v1 | 0.721 | 0.620 | 0.540 | 0.573 | 0.508 | 1.42 |
| FRA | 1.040 | 0.689 | 0.476 | 0.523 | 0.476 | 2.18 |
| GER (pinned) | 1.044 | 0.735 | 0.525 | 0.579 | 0.645 | 1.62 |
| USA (pinned) | 1.191 | 0.988 | 0.829 | 0.889 | 0.956 | 1.25 |
| KMV Table 3 | 1.153 | 0.888 | 0.697 | 0.754 | 0.70 | 1.65 |

The continuous closed form E_t[V(t)] = V_erg·2β/(ζ + 2β) per component gives
1.000 for PL v2 (chain: 0.982), so the discretization does not distort this.

**Readings.**
1. **Entry at z = 0 under perpetual youth reduces PL v2 from 1.38 to 0.98, not
   to 0.72.** The reason: with ζ = 1/180, (1 − ζ)^144 = **45%** of the
   population is older than 36 working years, and Poland's persistent
   component (half-life 27 years) keeps building variance for that long.
   Only truncating ages at 36 years (0.713) brings the cross-section to the
   target. Perpetual youth cannot do that without an age state variable,
   i.e. deterministic retirement/death after 144 quarters, which is a
   different model from KMV. The remaining gap between "uniform 0–36"
   (0.78) and the target (0.72) is time aggregation (quarterly vs annual
   earnings); 0.713 lands near 0.72 because the exponential age weights
   favour the young.
2. **PL v2's ergodic overshoot is not an outlier.** FRA has a *larger*
   ergodic/target ratio (2.18) and calibrates. GER 1.62 and KMV themselves
   1.65 also calibrate. What breaks PL v2 is the **level** of ergodic
   variance (1.38, the highest of all chains here, above USA's 1.19, which
   also fails), not its mismatch to the data. With z = 0 entry, PL v2 would
   reach 0.98, about the FRA/GER ergodic level (1.04) at which the one-asset
   model calibrates. So the change would likely restore calibration. That
   is a conjecture, not checked: it requires changing `death.py` and running
   `swap_death`.
3. **The change would apply to all countries, not just Poland.** It lowers
   cross-sectional dispersion for every chain (FRA 1.04 → 0.69, GER 1.04 →
   0.74). Their calibrations (§20, STEP2_LOG S5) would have to be redone, and
   the KMV baseline (S1, which reproduces KMV exactly) would cease to be
   comparable.

**Assessment.**
- **Consistency with KMV: no.** KMV explicitly choose ergodic draws.
  Entry at z = 0 is a deviation from the model being replicated, and a
  visible one: it changes the earnings distribution of newborns and so the
  hand-to-mouth fraction, which is exactly what the calibration targets.
- **Justification: partial.** For entry at z = 0: (i) the process is
  *estimated* on a panel starting from z = 0, so a z = 0 entry makes the HANK
  cross-section closer to the population on which it was estimated; (ii)
  "newborn = worker entering the market" is a natural reading of perpetual
  youth with a 45-year lifetime. Against: (i) it does not close the gap
  (0.98 vs 0.72) as long as death is Poissonian; the full fix requires an
  age state variable, which changes the model class; (ii) deviating from
  KMV for **one** reason (Poland does not calibrate) looks like
  target-shooting unless it is applied uniformly and justified ex ante;
  (iii) KMV themselves accept the same type of gap (1.07–1.15 vs 0.70).
- **Recommendation.** Keep the baseline consistent with KMV (ergodic draws).
  Treat z = 0 entry as a **robustness variant for all four countries**
  together, with this table as motivation, not as a fix for Poland. The
  primary decision about Poland's one-asset vs two-asset model (§24.6) should
  not depend on this variant. To discuss with the supervisor, together with
  the open question about two-asset.

---

## Open decisions (not yet settled)

- **One-asset vs two-asset for Poland — REOPENED 2026-10-20 (§24.6).** On the
  v2 chain (`var_log_earns` from `std_loginc`) the one-asset model + death
  does **not** calibrate: MPC max ≈ 0.07 vs 0.20, HtM 6%. Robust to 9 starting
  points. The §20 "RESOLVED" (below) applies to the v1 targets. Before
  deciding on two-asset (STEP2_LOG §S6, ≈ 2 weeks): (i) check whether a chain
  from the v2 profile floor with lower var₂ (`output/pol_v2/var2_1.0/`, inside
  the noise band) calibrates; if so, there is a defensible one-asset route with
  var₂ chosen within the statistical uncertainty, which would need a
  justification for the supervisor; (ii) a question for the supervisor: is
  something like the GRID-USA treatment ("documented sensitivity") acceptable
  for Poland? Poland is the main country of the thesis, so probably not — to
  be confirmed.
- **Robustness test for the level target (§24.3/§24.4): `std_resid` + model
  moments on ages 25–55.** The baseline (§24.4) uses raw `std_log_inc²` and a
  36-year panel; the two biases offset to within a few percent. The fully
  consistent alternative is GRID `std_resid²` (`data/std resid.csv`, already
  available for all four countries) together with model moments computed only
  on panel columns for ages 25–55 (31 years). It requires re-estimating all four
  countries. To discuss with the supervisor: is it worth it as a robustness test
  in the thesis, or is the §24.3 table enough. Before running it, also confirm
  GER's `nobs` inconsistency between exports (24.3).
- ~~**If `β₁` again runs to its bound after widening to 8.0:** fix it and
  estimate the remaining five.~~ **Resolved** (§8): `β₁` settled interior at
  4.27 (FRA) / 5.34 (GER) / 2.75 (USA), objective sweeps confirm. No action.
- **The `var Δ1y` gap** (§8d): implement the fix (discretize component 1 as a
  near-iid quarterly shock) at step S4. Frequency question is **settled** —
  step-2 model is quarterly (STEP2_LOG §S2), our chain is quarterly, no `P⁴`.
- **`beta` calibration convention for the step-2 comparison** (STEP2_LOG open
  list): re-calibrate `beta` per chain, or hold it at the baseline value?
  Re-calibration is the cleaner validation design. Supervisor question.
- **GER `β₂` calibration** (§8b). (i) ~~stationary vs lifecycle~~ **settled**:
  the step-2 code is infinitely-lived / ergodic, so the restriction is needed.
  Still open: (ii) whether the ergodic-sd anchor should be Germany's own observed
  dispersion rather than KMV's US value; (iii) whether to add a ninth identifying
  moment (longer-horizon autocovariance) as a robustness check. Both for the
  supervisor.
- **USA `β₂` pinned — DONE (2026-08-30), but it did not fix step 2.** USA's
  unrestricted estimate implied ergodic `var(log e) ≈ 1.79` vs KMV's own 1.07,
  so USA was re-estimated with `β₂` derived from the same `ErgodicVarConstraint`
  as GER (→ `β₂ = 0.0136`, `σ₂` 1.32 → 1.53, ergodic sd of the persistent
  component = KMV's 0.954). DE and TikTak agree; OOS f 0.033 (vs 0.006
  unrestricted — the constraint costs USA fit, as it does GER). `output/usa/`
  now holds the pinned chain; the unrestricted one is kept at
  `output/usa_free_beta2/`. **This is now the USA canonical.** But in step 2
  USA *still* does not calibrate in one-asset even with `β₂` pinned + stochastic
  death (STEP2_LOG §S5): the GRID-USA chain lacks KMV's leptokurtic tail
  (`kurt Δ1y` 12.8 vs 17.8), so too few households reach the constraint whatever
  the discount factor. FRA and GER *do* calibrate with death + β-het. Proposed:
  the S1 baseline (KMV's chain) is the US validation; GRID-USA is a sensitivity.
- **The §8b ergodic-plausibility argument needs a light revisit** given §S5:
  it is still right that a stationary HANC confronts households with the ergodic
  distribution, but the reference point is KMV's 1.07, not "implausible for any
  labour market". GER's pinned ergodic sd 0.95 (var 0.91 for the persistent
  component) is right at KMV's; the argument is "match KMV's ergodic dispersion
  because the German panel can't identify `β₂`", which is cleaner than the
  original "avoid implausible values" framing. Reword when convenient.
- **Restoring a genuine second optimiser (§9) — optional, not blocking.** The
  line scan + DE cross-run stability already carry the identification argument,
  so the estimates are settled. But if a stronger form of the evidence is wanted,
  swapping TikTak's Nelder-Mead local phase for a pattern-search or Powell method
  (or running its objective on a larger, smoother simulation) would plausibly fix
  the stalling and give back a real independent check. Supervisor's call whether
  it is worth doing.
- **Whether the German near-permanent-shock result** (weaker `β₂` decay than
  USA/FRA in the unconstrained fit — see §8b) **has an institutional
  interpretation** worth developing in the text (sectoral bargaining, employment
  protection), or should be reported without one.
- **Chain size for step 2:** 75 states is the current default; may need to drop
  to 33 if the model solve is too slow.
- ~~**One-asset vs two-asset model in step 2 — Poland is the trigger.**~~
  **[REOPENED 2026-10-20 — see the first item in this list and §24.6; the
  'resolved' below holds only for the v1 std_perm targets.]**
  **RESOLVED 2026-10-12, see §20: Poland calibrates cleanly, one-asset + death
  is sufficient.** The rule-of-thumb threshold below (`β₁ ≲ 3`) would have
  pointed the wrong way for Poland (`β₁ = 2.85`) — the actual `swap_death.py`
  run is decisive and settles it instead. Two-asset remains wanted on other
  grounds (cross-sectional MPC distribution, wealthy-HtM mechanism, portfolio
  channel) if the thesis needs those, just not *forced* by a calibration
  failure. Original rule-of-thumb text kept below for the record:
  the one-asset ARS model + KMV stochastic death (`death.py`,
  STEP2_LOG §S5) calibrates cleanly for **FRA and GER** — MPC 0.20, monetary IRF
  within 7% of baseline, KMV's 20/80 split reproduced (§S3). **USA** does not
  calibrate (thin GRID-USA tail, §12) but is validated by the S1 baseline on
  KMV's chain. So one-asset is the plan.
  The trigger for switching to two-asset was going to be: the Polish estimated
  chain is GRID-USA-shaped — `kurt Δ1y` well below ≈ 15 *and/or* a slow-decaying
  transitory component (`β₁` ≲ 3). If Poland looks like France/Germany
  (`kurt Δ1y` ≳ 15, `β₁` ≳ 4), one-asset + death handles it. If it looks like
  GRID-USA (thin tail, persistent transitory), the one-asset model will not
  calibrate for Poland and the two-asset model (STEP2_LOG §S6, ≈ 2 weeks, SSJ
  `examples/two_asset`) is required.
- **Aggregate match vs channel-level match (STEP2_LOG §S3).** The 20/80
  direct–indirect *split* matches KMV, but the *composition* of the indirect
  effects does not (ARS asset-return/capital-gains channel +22% vs KMV's −2%;
  transfer channel 16% vs 32%) — because of ARS's `capitalization` block and
  fast fiscal rule vs KMV's φ = ω̄ profit-neutralisation. If the Polish results
  report a decomposition *by channel* rather than just direct/indirect, that
  comparison to KMV will not be clean without modifying the non-household
  blocks. Is the aggregate match sufficient for the thesis? Supervisor question.
- **Polish `beta` calibration target.** USA/FRA/GER calibrate `beta` to known
  asset targets (A = 20 = 500% of annual GDP; SCF Lorenz). Poland's wealth
  moments (NBP HFS / ECB HFCS) are not in hand, and Polish household wealth
  differs from the US in level *and* composition — so Polish `beta` lands
  elsewhere and that alone moves MPCs and transmission. "What do we calibrate
  Polish `beta` to, and what is the asset target" is a thesis-design question
  for Zoch, not an implementation detail. Raise before the Polish stage.
- ~~**Poland:** earnings targets (Polish panel data / EU-SILC) and wealth
  moments (NBP HFS / ECB HFCS) still pending.~~ **Resolved** — Polish GRID
  targets built for 2000-2020 (§18), wealth moments still pending and still a
  Żoch question.
- ~~**Is the canonical GER point the floor of the valley?**~~ **Resolved
  2026-10-12, see §9's extended-line-scan addendum:** yes — the objective
  rises sharply and monotonically past the DE point for both GER and FRA
  (GER quadruples, FRA rises over 40x by 100% further along the same
  direction). The canonical point is a genuine local minimum, not a point on
  a still-descending slope.
