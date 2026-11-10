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

**Erratum (2026-11-09, §29.7).** This verification does not hold. With
correct units every reset implementation gives var Δ1 ≈ 0.23 at Table 3
parameters. The 0.40 is reproduced only with rates 4× too fast, most likely
a units slip in the original check (not kept). The jump type matters for the
level and Δ5 through the sampling design (lifecycle vs ergodic), not for
var Δ1. KMV's own code (replication package) uses resets. The additive
choice stands because the working paper (NBER w21897 §4.3, CEP 2016 eq. 31)
writes the process additively. The reset version is a robustness check
(§29.8).

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

### 24.12 Summary of §24 (status as of 2026-10-20)

**Done.**
- Problem diagnosed (24.1). `std_perm` = log of the 3-year average,
  age-residualized, on sample [PE]. Target corrected to `std_loginc²` = 0.722
  (24.2). Decision to keep raw `std_log_inc` + the 36-year panel for all
  countries (24.3–24.4).
- Poland v2 re-estimated (DE + TikTak), var₂ profile, report (24.5).
  `swap_death` on the v2 chain: **does not calibrate** (24.6).
- Re-runs on v2: sub-periods (24.8), nested specs (24.9), recovery (24.10),
  bootstrap (24.11).
- §25: newborns in HANK (KMV = ergodic draw; z = 0 would give 0.98, not
  0.72). `death.py` unchanged.

**Do any of these results change the §24.5 conclusions? No. They strengthen
them.**
1. *Does Poland have the lowest dispersion?* Still **no, the highest** among
   the canonical estimates. The bootstrap CI for the persistent ergodic sd is
   [1.00, 1.13], above FRA 0.90 and the pinned USA/GER 0.95. Both sub-periods
   give it even higher (1.11 / 1.17).
2. *Do the data prefer a value below 0.91?* Still **no, above it**. The var₂
   bootstrap CI [1.00, 1.28] excludes 0.91 (stronger than the interpolation
   in 24.5).
3. *Is there a clear floor?* Still **yes**. Recovery: β₂ point-identified
   (−2%, 4/4 reps, no v1-style bimodality). Nested V3 (β₂ = 0) rejected ×3.1.
   The bootstrap CI lies inside the profile's noise band. The new nuance is
   weak precision on the right side (asymmetric β₂ CI, sub-periods drifting
   toward higher var₂).
- One genuinely new observation: in v2 the persistent component is less
  stable across sub-periods than in v1 (half-life 29 vs 40 years pre/post).
  It is consistent with weak precision on the right side, not with a break.

**Left (not run, waiting on a decision).**
- **One-asset vs two-asset for Poland** (Open decisions, 24.6). Inputs:
  `swap_death` on `output/pol_v2/var2_1.0/` (not run, as instructed), the
  §25 variant (z = 0 entry, all countries), the two-asset model.
- Exact f at var₂ = 0.9103 on the v2 targets (one run, ~35 min). Now low
  priority, since the bootstrap CI already excludes 0.91.
- Robustness test with `std_resid` + ages 25–55 (Open decisions, 24.3/24.4)
  and confirming GER's `nobs` inconsistency between exports.
- `swap_death` on the v2 sub-periods (deliberately not run, 24.8).
- Updating the supervisor note (tag `notatka-promotor-v1`): the Polish
  §18–§23 results are to be marked as v1, superseded by §24.

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

## 26. Poland v3: residualized level, ages 25–55, trends in the moments (2026-11-08)

**Scope after the meeting with the supervisor.** From here on only Poland is
re-estimated. FRA/GER/USA stay exactly as they are (method validation, §0a/§11)
and are **not** recomputed under the v3 definitions. Nothing from v1/v2 is
overwritten: v3 writes `targets/pol_grid_v3*.json`, `output/pol_v3/`,
`output/pol_v3_end/` and new scripts. For Poland only, this implements "option
1" from §24.3 (the robustness item in Open decisions).

### 26.1 Targets v3 (A.1)

`targets/pol_grid_v3.json` = v2 with
`var_log_earns = mean_{2000..2020}(std_resid²) = 0.69119` (GRID `std_resid`,
`data/std resid.csv`, POL, Male, 25–55; age-residualized annual log earnings
on sample [M], `nobs_resid = nobs_loginc` in every year, checked again in the
script). The other seven moments are copied from v2 unchanged (they are
residualized changes already). Script: `scripts/pol_v3_targets.py`.

One small inconsistency found while rebuilding the year panel: the three
fractions in the targets are interpolated from the **across-year mean
percentiles** (`grid_loader.fractions_from_percentiles`), whereas the year panel
(`year_panel_pol_v2.csv`, used by the bootstrap since §19) interpolates **year by
year**, and the mean of those differs by +0.0059 / −0.0013 / +0.0015
(<10 / <20 / <50%). Variances and kurtoses match to machine precision. Decision:
targets stay as in v2 (the brief says so, and it keeps v2 vs v3 a one-moment
comparison). The gap is ≤ 1.3% of the target and is noted wherever the year
panel is used (trends, bootstrap, Ω_data).

### 26.2 Model moments on ages 25–55 only (A.2)

`moments_from_panel(logy, n_age_cols=None)` and
`model_moments(params, n_age_cols=None, **sim_kwargs)`. With `n_age_cols=31`
only the first 31 columns of the 36-year panel (ages 25–55) are used, for the
level and for 1- and 5-year changes (5-year changes therefore only pair ages
25→30 … 50→55). The default `None` is **bit-for-bit unchanged**:
`tests/test_age_window.py` compares against a verbatim copy of the old function.
A separate check against the `HEAD` version of `simulate.py` gave identical
dicts at KMV Table 3 and at the PL v2 estimate. The first k simulated years do
not depend on how many years follow (same RNG calls in the same order), so
`model_moments` simulates only 31 years when `n_age_cols=31`. The test confirms
that the moments equal the 36-year panel sliced to 31 columns exactly (≈14% less
work per evaluation). The option is threaded through `discretized_moments`,
`discretize_process(match_var_log=True)` and `run.build_table`. All Poland v3
runs use `n_age_cols=31`.

### 26.3 Linear trends in the eight moments (A.3)

Annual panel `output/pol_v3/year_panel_pol_v3.csv` = `year_panel_pol_v2.csv`
with `var_log = std_resid²` (v2's `std_loginc²` kept as `var_log_v2`).
OLS of each series on a constant and the year. Newey–West slope SE, Bartlett
kernel, lag 2 (= ⌊4(T/100)^{2/9}⌋ for T = 16–21; no small-sample correction).
Plain OLS SEs are in `trends.csv`.

**How the different series lengths are treated.** GRID indexes a change by its
base year: Δ1 in year t is t→t+1 (2000–2019), Δ5 in year t is t→t+5
(2000–2015). The level runs 2000–2020. Each moment's trend is fitted on **its
own available years**, with no padding and no common sub-window. Its "window
midpoint" is the mean of its own years (2010 / 2009.5 / 2007.5). Checked in the
script: the trend at that midpoint equals the series mean to machine
precision (an OLS identity with a constant). So the v3 baseline targets (window
means) are exactly "the trend at the midpoint"; for the fractions this holds up
to the §26.1 interpolation gap.

| moment | years | mean (= trend at mid) | slope / yr | NW se | NW t | trend 2000 | mid year | trend at series end | trend 2019 |
|---|---|---|---|---|---|---|---|---|---|
| var log earns (`std_resid²`) | 2000–2020 | 0.6912 | −0.00671 | 0.00168 | −4.0 | 0.758 | 2010 | 0.624 (2020) | 0.631 |
| var Δ1 | 2000–2019 | 0.2595 | −0.00385 | 0.00060 | −6.4 | 0.296 | 2009.5 | 0.223 (2019) | 0.223 |
| var Δ5 | 2000–2015 | 0.4765 | −0.01026 | 0.00106 | −9.6 | 0.554 | 2007.5 | 0.400 (2015) | 0.359* |
| kurt Δ1 | 2000–2019 | 15.24 | +0.098 | 0.054 | +1.8 | 14.31 | 2009.5 | 16.17 (2019) | 16.17 |
| kurt Δ5 | 2000–2015 | 9.878 | +0.150 | 0.028 | +5.3 | 8.75 | 2007.5 | 11.01 (2015) | 11.61* |
| frac Δ1 < 10% | 2000–2019 | 0.4558 | +0.00245 | 0.00196 | +1.2 | 0.433 | 2009.5 | 0.479 (2019) | 0.479 |
| frac Δ1 < 20% | 2000–2019 | 0.6681 | +0.00269 | 0.00101 | +2.7 | 0.642 | 2009.5 | 0.694 (2019) | 0.694 |
| frac Δ1 < 50% | 2000–2019 | 0.8445 | +0.00138 | 0.00045 | +3.1 | 0.831 | 2009.5 | 0.858 (2019) | 0.858 |

\* linear extrapolation 4 years past the last data year (2015).

Files: `output/pol_v3/trends/{trends.csv, trends.json, trend_residuals.csv,
trends.png, v3_end_meta.json, log.txt}`.

Readings:
- **Dispersion falls over the window, and tails get fatter.** All three
  variances fall significantly (−19% / −28% / −32% of the mean over the window).
  The kurtosis of 5-year changes rises (+23%), and there are more small changes.
  Polish earnings risk became less dispersed but more concentrated in rare,
  large changes.
- **A linear trend is a rough description.** In the kurtoses and fractions the
  plot shows a level shift around 2008–2009 rather than a steady drift (and a
  reversal after 2016 in Δ1). NW t-statistics with T = 16–21 and lag 2 are also
  only indicative. This matters for v3_end (below) and for the detrended
  bootstrap (B.4 ii), whose residuals partly carry the break.

### 26.4 v3_end targets: the trend at the end of the window (A.4, robustness only)

`targets/pol_grid_v3_end.json` = trend value at base year **2019** (last year
with 1-year changes) for every moment, written as
`v3 target + slope × (2019 − own midpoint)`. For the five non-fraction moments
this is exactly the trend line. For the fractions it keeps the target's
mean-percentile construction (§26.1) and adds the trend shift, so `v3_end − v3`
is the trend alone.

- **Extrapolated:** var Δ5 and kurt Δ5 (data end in base year 2015), 4 years
  past the data. Marked in `output/pol_v3/trends/v3_end_meta.json`.
- **Not extrapolated:** the level (data to 2020) and the four Δ1 moments.
- Alignment caveat: by the year a change *ends*, the last Δ5 (2015→2020) and
  the last Δ1 (2019→2020) end in the same year. By the change's midpoint they
  are 2 years apart (2017.5 vs 2019.5). Following the brief, the base-year
  convention is used; it is the most conservative reading (4 years of
  extrapolation).

**Targets v2 vs v3 vs v3_end**

| moment | v2 | v3 | v3_end | v3_end vs v3 |
|---|---|---|---|---|
| var log earns | 0.7216 (`std_loginc²`) | **0.6912** (`std_resid²`) | 0.6308 | −8.7% |
| var Δ1 | 0.2595 | 0.2595 | 0.2229 | −14.1% |
| var Δ5 | 0.4765 | 0.4765 | 0.3585* | −24.8% |
| kurt Δ1 | 15.24 | 15.24 | 16.17 | +6.1% |
| kurt Δ5 | 9.878 | 9.878 | 11.61* | +17.5% |
| frac Δ1 < 10% | 0.4499 | 0.4499 | 0.4732 | +5.2% |
| frac Δ1 < 20% | 0.6694 | 0.6694 | 0.6950 | +3.8% |
| frac Δ1 < 50% | 0.8431 | 0.8431 | 0.8562 | +1.6% |

v3 differs from v2 in one target (−4.2%) and on the model side (ages 25–55).
v3_end describes a different period altogether; it is a robustness test, not
the baseline.

### 26.5 Estimation of Poland v3 (B.1)

Script `scripts/pol_v3_estimate.py`. Setup:
- estimation panel: 50k workers, ages 25–55 (`n_age_cols=31`), CRN;
- DE: 5 seeds (7, 101, 202, 303, 404), maxiter 200, popsize 15, NM polish 800;
- TikTak: n_sobol 4096, n_local 30;
- warm start at KMV Table 3;
- every start re-evaluated out of sample on a 100k panel (seed 987654).

The canonical estimate is the start with the lowest OOS f. Per-start results
are in `output/pol_v3/runs/*.json`, and the summary is `output/pol_v3/report.json`.

**Spread across the six starts** (rates quarterly; var₂ = λ₂σ₂²/(2β₂); hl₂ =
half-life of z₂ in years):

| start | f in-sample | f OOS | λ₁ | β₁ | σ₁ | λ₂ | β₂ | σ₂ | var₂ | hl₂ |
|---|---|---|---|---|---|---|---|---|---|---|
| **DE seed 404** | 0.00734 | **0.00762** | 0.1623 | 2.801 | 1.998 | 0.01004 | 0.005665 | 1.193 | 1.260 | 30.6 |
| DE seed 202 | 0.00734 | 0.00762 | 0.1623 | 2.802 | 1.999 | 0.01004 | 0.005667 | 1.192 | 1.259 | 30.6 |
| DE seed 101 | 0.00734 | 0.00762 | 0.1623 | 2.800 | 1.997 | 0.01004 | 0.005703 | 1.194 | 1.255 | 30.4 |
| DE seed 7 | 0.00735 | 0.00764 | 0.1624 | 2.801 | 1.998 | 0.01004 | 0.005665 | 1.193 | 1.260 | 30.6 |
| DE seed 303 | 0.00735 | 0.00764 | 0.1624 | 2.801 | 1.998 | 0.01004 | 0.005664 | 1.193 | 1.260 | 30.6 |
| TikTak | 0.00751 | 0.00773 | 0.1609 | 2.829 | 2.011 | 0.01062 | 0.005216 | 1.143 | 1.331 | 33.2 |

- **One basin.** The five DE seeds agree to 3–4 digits; their OOS f differ by
  < 0.3%.
- **TikTak lands on the same ridge, slightly worse.** It is 1.5% higher in OOS
  f and moved along the known β₂/σ₂ direction (β₂ −8%, var₂ +5.6%, hl₂ +2.6
  yr). This is the same pattern as §9/§18, the weakly identified persistent
  component.
- **No line scan was run.** The brief asks for one only if starts land in
  different basins. The best two (DE 404 and DE 202) differ by at most 0.07% in
  any parameter, so the trigger (> 10%) did not fire. The largest spread across
  all six starts is β₂, 7.9% (TikTak vs the best).
- **Not at a bound.** Log positions inside the bounds run from 0.50 to 0.80
  (λ₂ 0.50, β₂ 0.53).

**Fit of the eight moments** (OOS panel, canonical; discretized = 5×15 chain
with `match_var_log`):

| moment | target v3 | model | rel. error | chain 5×15 |
|---|---|---|---|---|
| var log earns | 0.6912 | 0.690 | −0.1% | 0.693 |
| var Δ1 | 0.2595 | 0.261 | +0.8% | 0.176 |
| var Δ5 | 0.4765 | 0.475 | −0.2% | 0.409 |
| kurt Δ1 | 15.24 | 14.51 | −5.4% | 17.02 |
| kurt Δ5 | 9.878 | 9.954 | +1.0% | 10.55 |
| frac < 10% | 0.4499 | 0.480 | +6.4% | 0.659 |
| frac < 20% | 0.6694 | 0.657 | −2.0% | 0.676 |
| frac < 50% | 0.8431 | 0.837 | −0.9% | 0.831 |

The misfit pattern is the same as in v2: kurt Δ1 too low, too many small
changes. The chain's known biases (§8d: var Δ1 −32%, mass at < 10%) are
unchanged.

**v2 vs v3 (canonical estimates).** v2 is the DE estimate from
`output/pol_v2/params.json`, on a 36-year panel.

| | v2 | v3 | change |
|---|---|---|---|
| λ₁ | 0.1623 | 0.1623 | 0.0% |
| β₁ | 2.832 | 2.801 | −1.1% |
| σ₁ | 2.001 | 1.998 | −0.2% |
| λ₂ | 0.00987 | 0.01004 | +1.7% |
| β₂ | 0.006376 | 0.005665 | −11.2% |
| σ₂ | 1.209 | 1.193 | −1.3% |
| OOS f | 0.00728 | 0.00762 | (different targets and age window: not comparable) |
| var₁ (ergodic, additive formula) | 0.115 | 0.116 | |
| **var₂** | **1.132** | **1.260** | **+11.3%** |
| ergodic sd z₂ | 1.064 | 1.122 | +5.5% |
| ergodic sd z₁ + z₂ | 1.117 | 1.173 | |
| half-life z₁ (yr) | 0.061 | 0.062 | |
| half-life z₂ (yr) | 27.2 | 30.6 | +3.4 yr |
| years between persistent shocks | 25.3 | 24.9 | |
| chain ergodic var(log e), 5×15 exported | 1.38 | 1.51 | +9.6% |

Chain values are computed from the exported `income_process_{grid,pi}.txt`.
`match_var_log` scales the persistent component's chain variance above the
formula value (1.386 vs 1.260 for z₂ in v3, §8c).

**What changed relative to v2, and why.**
- The transitory component is unchanged. Its moments (var Δ1, kurt Δ1,
  fractions) are the same targets in v2 and v3.
- The level target fell 4.2% (0.722 → 0.691), yet **var₂ rose 11%**. The
  model side changed too: v3 averages the cross-sectional variance over ages
  25–55 only, i.e. without the five oldest, most dispersed ages of a panel that
  starts from z = 0. The average variance per unit of var₂ is therefore lower,
  and matching a level only 4% lower needs a larger ergodic variance and a
  slower drift (β₂ −11%).
- So v3 does not soften the v2 conclusion. Poland's long-run dispersion is
  higher in v3: var₂ 1.26, and the exported chain has ergodic var(log e) 1.51
  (persistent component 1.39).
- **KMV reference (changed 2026-11-09, §29.9).** The reference is the
  ergodic distribution of the **published KMV chain**
  (`annual-review/inputs/kmv_process`, the chain inside their model and our
  S1 HANK baseline): var(log e) **0.846** (sd 0.92), persistent component
  **0.647**.
  - Earlier text used 0.91, our additive formula applied to KMV's parameters,
    and 0.655, the reset formula at the rounded Table 3 parameters. Both are
    superseded.
  - Against the chain, Poland v3 is 1.51 vs 0.85 in total and 1.39 vs 0.65 for
    the persistent component.

### 26.6 Bootstrap of Poland v3 (B.4)

Script: `scripts/pol_v3_bootstrap.py`. Output: `output/pol_v3/bootstrap/`.

**Setup.** 200 reps per design, with the same reduced budget as §19/§24.11:
- 5k workers, DE 10×5, NM 50;
- warm start at the canonical estimate;
- ages 25–55.

**Two designs:**
- **(i) `years`**: resample the 21 years with replacement (§19 design). Each
  moment is the mean over the years drawn that have that moment, so the
  shorter Δ5/Δ1 series carry fewer draws.
- **(ii) `detrended`**: resample whole year vectors of the detrended residuals
  (`trend_residuals.csv`, §26.3) and add them to the trend at its midpoint
  (= the v3 target). This removes the 2000–2020 trend from the variability, so
  only fluctuations around the trend remain.

**Validation of the reduced budget.** 3 resamples per design, re-estimated at
both budgets (reduced vs §19's reference budget: 15k workers, DE 20×6, NM 100), in
`validation.json`. Deviations, reduced vs reference:

| | relative deviation |
|---|---|
| σ₁, σ₂ | ≤ 1% |
| λ₁, β₁, λ₂ | −1% to −4%, mean −2% |
| **β₂** | **0% to −12%, mean −5.6%** |

Every deviation is below §19's 15% threshold. **But the β₂ deviation is
systematic.** The reduced budget lands lower on the β₂/σ₂ ridge, i.e. at a
higher var₂. This shows up as an upward shift of the bootstrap means (var₂
1.32–1.33 vs point estimate 1.26), and §24.11 had the same shift (β₂ −3.7%).
Below I therefore report the SEs and **recentred** percentile intervals
(percentiles shifted by mean − point). The raw percentiles are in
`bootstrap_summary.csv`.

**Resampled targets.** The SD across reps is smaller under (ii) than under (i):

| moment | (i) years | (ii) detrended |
|---|---|---|
| level | 0.0102 | 0.0062 |
| var Δ5 | 0.0126 | 0.0038 |
| kurt Δ5 | 0.204 | 0.105 |
| var Δ1 | 0.0058 | 0.0026 |
| kurt Δ1 | 0.253 | 0.211 |

So most of the year-to-year variability in the variances is the trend itself.
The kurtosis of Δ1 is mostly noise around the trend.

**Results**

| | point | SE (i) years | rel. | SE (ii) detrended | rel. |
|---|---|---|---|---|---|
| λ₁ | 0.1623 | 0.0038 | 2.4% | 0.0026 | 1.6% |
| β₁ | 2.801 | 0.066 | 2.4% | 0.059 | 2.1% |
| σ₁ | 1.998 | 0.011 | 0.5% | 0.0073 | 0.4% |
| λ₂ | 0.01004 | 0.00042 | 4.2% | 0.00025 | 2.5% |
| β₂ | 0.005665 | 0.00063 | **11.1%** | 0.00031 | **5.6%** |
| σ₂ | 1.193 | 0.019 | 1.6% | 0.010 | 0.9% |
| var₂ | 1.260 | 0.097 | 7.7% | 0.050 | 3.9% |
| ergodic sd z₂ | 1.122 | 0.042 | 3.7% | 0.022 | 1.9% |
| half-life z₂ (yr) | 30.6 | 3.7 | 12.0% | 1.9 | 6.1% |

**95% intervals for the persistent component** (recentred percentiles; raw
ones in brackets):

| | (i) years | (ii) detrended |
|---|---|---|
| var₂ | **[1.06, 1.44]** ([1.12, 1.50]) | **[1.17, 1.35]** ([1.24, 1.41]) |
| ergodic sd z₂ | [1.03, 1.20] | [1.08, 1.16] |
| half-life z₂ (yr) | [23.1, 36.8] | [26.9, 34.1] |

**Readings.**
- **(ii) halves the SEs of the persistent component compared with (i).** Under
  (i), resampling years treats the 2000–2020 decline in dispersion as noise. A
  draw heavy in early years has a high level and a high var Δ5. Under (ii) only
  deviations from the trend remain. Which one is "right" depends on what the
  target is meant to be:
  - if it is the mean of a stationary process observed over 21 years, use (i);
  - if it is the midpoint of a trending process, use (ii).
  - The trend is strongly significant (§26.3), so (ii) is the more coherent
    description of uncertainty around **this** target. (i) is the conservative
    upper bound. Both are reported.
- **The two designs agree on the transitory component** (SEs of 0.4–2.4%).
  They differ only in how they weigh the persistent one.
- **Both intervals lie far above KMV's persistent component** (0.65 in the
  published KMV chain, §29.9). The lower ends are 1.06 and 1.17, and v2's CI
  was [1.00, 1.28]. They also exclude the superseded anchor 0.91.
  - Caveat: Poland's var₂ is a formula value and KMV's 0.65 a chain value.
  - The chain sits about 10% above the formula for Poland (1.39 vs 1.26) and
    7% above for KMV (0.647 vs the exact reset formula 0.603), so the
    comparison is roughly like-for-like in scale.
- Caveats, as in §19 and §24.11:
  - this is year-to-year variability of the GRID statistics, not microdata
    sampling error;
  - every rep is warm-started at the point estimate, so the intervals are
    conditional on its basin. §26.5 found one basin.
  - The years are treated as independent, but the residuals carry the
    2008–09 level shift, i.e. they are serially correlated. A block bootstrap
    would widen (ii). It was not run; T = 16–21 years is too short for
    meaningful blocks.

### 26.7 Profile, recovery, nested specs for Poland v3 (B.3)

Script: `scripts/pol_v3_robustness.py`. Outputs:
- `output/pol_v3/{var2_profile, recovery, nested_specs}/`;
- `output/pol_v3/robustness_summary.json`.

Budget as in §21–§23: DE 100 × 12, NM 400, 50k workers, ages 25–55, OOS on
the 100k panel (seed 987654).

**var₂ profile.** var₂ = λ₂σ₂²/(2β₂) is pinned at m × 1.260 and the other
five parameters are re-estimated, starting from the canonical estimate. A
5 × 15 chain is exported for every point (`var2_profile/chain_*`), ready for
step 2.

| m | var₂ | ergodic sd z₂ | OOS f | f / f_min | worst moments |
|---|---|---|---|---|---|
| 0.6 | 0.756 | 0.869 | 0.02178 | ×2.86 | var log −10.8%, var Δ5 +6.8% |
| 0.8 | 1.008 | 1.004 | 0.00906 | ×1.19 | frac<10 +6.5%, kurt Δ1 −5.2%, var Δ5 +3.0% |
| 0.9 | 1.134 | 1.065 | 0.00787 | ×1.03 | (as at the optimum) |
| **1.0 (control)** | **1.260** | **1.122** | **0.00762** | 1 | frac<10 +6.4%, kurt Δ1 −5.4% |
| 1.1 | 1.386 | 1.177 | 0.00776 | ×1.02 | (as at the optimum) |
| 1.25 | 1.575 | 1.255 | 0.00829 | ×1.09 | var Δ5 −2.4% |
| 1.5 | 1.890 | 1.375 | 0.00921 | ×1.21 | var Δ5 −3.8% |

- **The control point reproduces the canonical f exactly** (0.00762), so the
  pinned re-estimation finds the optimum again.
- **The profile has a clear minimum at 1.0×.** It is flat within ±10%
  (+2–3% in f), and the level and var Δ5 start to bind outside it.
- **It is asymmetric.** It rises steeply toward low var₂: at 0.8× f is +19%,
  at 0.6× ×2.9, with the level 11% too low. Toward high var₂ it is flatter:
  +9% at 1.25×, +21% at 1.5×.
- **The KMV reference lies outside the profiled range.** The persistent
  component of the published KMV chain (0.65, §29.9) is ≈ 0.51× Poland's
  var₂. That is below the lowest profile point (0.6×, already f ×2.9, level
  −11%).
  - The superseded anchor 0.91 (≈ 0.72×) would sit on the steep side at
    ×1.5–1.8 the minimum.
- Compared with v2 (§24.5), the shape is the same, shifted up by roughly the
  11% change in the point estimate.

**Recovery.**
- Design: "true" parameters = canonical v3. The targets are the model moments
  on a panel with a different seed (101 or 202) and size (50k or 300k).
  Re-estimation starts from **KMV Table 3**, i.e. far from the truth.

| rep | λ₁ | β₁ | σ₁ | λ₂ | β₂ | σ₂ | ergodic sd z₂ (true 1.122) | f |
|---|---|---|---|---|---|---|---|---|
| 101, 50k | +4.5% | +7.1% | −0.4% | +0.9% | −0.2% | −0.8% | 1.120 | 2e-5 |
| 101, 300k | +0.9% | +0.5% | 0.0% | −1.2% | −1.5% | +0.3% | 1.128 | 0 |
| 202, 50k | +1.3% | +0.1% | −0.4% | −3.7% | +0.3% | +1.9% | 1.121 | 0 |
| 202, 300k | +4.6% | +6.0% | −0.4% | −0.7% | −1.5% | +0.4% | 1.132 | 2e-5 |

- **4 / 4 recover the true basin.** β₂ is within ±1.5% and the ergodic sd of
  z₂ within +1%. The same holds in v2 (4/4); v1 had 1 failure in 4.
- **The persistent component is recovered better than the transitory one.**
  λ₁ and β₁ miss by up to 4.5% and 7%, which is the size of their asymptotic
  SE (§27, 5–6%). This supports §27's view that the bootstrap understates the
  uncertainty of (λ₁, β₁): those parameters really do move by 5% with
  simulation noise alone.

**Nested specifications** (re-estimated with the same budget):

| variant | free params | OOS f | × full | worst moment | v2: × full |
|---|---|---|---|---|---|
| full | 6 | 0.00762 | 1 | frac<10 +6.4% | 1 |
| V1 single component | 3 | 0.460 | ×60 | var log −59%, frac<10 +25% | ×66 |
| V2 deterministic arrival (AR(1)) | 2 | 2.34 | ×307 | kurt Δ1 −80%, kurt Δ5 −70% | ×324 |
| V3 β₂ = 0 (random walk) | 5 | 0.0166 | **×2.2** | var Δ5 −8.8%, var log +3.5% | ×3.1 |

- **The conclusions are as in v2 (§24.9).**
  - Without jumps (V2) the kurtoses cannot be fitted at all.
  - With one component (V1) the level and the 1-yr fractions cannot both be
    fitted.
  - Of the restrictions, only the random walk (V3) is close to the full model.
- V3 is relatively better in v3 than in v2 (×2.2 vs ×3.1): with the window cut
  at 55, the persistent component's mean reversion matters less. It is still
  rejected. f is 2.2× higher, var Δ5 is 9% too low and the level 3.5% too
  high, the profile signature of β₂ → 0.

### 26.8 v3_end: estimation on the end-of-window trend targets (B.2)

Setup: `scripts/pol_v3_estimate.py --version v3_end`, 3 DE seeds
(7, 101, 202) plus TikTak, same budget as B.1. Outputs:
`output/pol_v3_end/{runs/, report.json, params.json, table2.csv}`.

| start | f in | f OOS | λ₁ | β₁ | σ₁ | λ₂ | β₂ | σ₂ |
|---|---|---|---|---|---|---|---|---|
| **TikTak** | 0.00663 | **0.00691** | 0.1577 | 2.821 | 1.950 | 0.00679 | 0.000104 | 1.126 |
| DE 101 | 0.00654 | 0.00694 | 0.1601 | 2.900 | 1.949 | 0.00683 | **0.0001** | 1.124 |
| DE 202 | 0.00654 | 0.00698 | 0.1624 | 2.985 | 1.950 | 0.00683 | **0.0001** | 1.124 |
| DE 7 | 0.00654 | 0.00699 | 0.1625 | 2.985 | 1.950 | 0.00683 | **0.0001** | 1.124 |

**β₂ goes to its lower bound (10⁻⁴) in every start.** The three DE seeds
sit exactly on it, and TikTak 4% above it. The persistent component becomes
a **random walk with jumps** (half-life about 1700 years). The formula value
var₂ = λσ²/(2β) ≈ 41 is meaningless: there is no ergodic distribution.

- **No line scan was needed.** All four starts are in the same corner. The
  best two differ by at most 4% (β₁).
- **No chain could be exported.** `discretize_process` correctly refuses a
  component with no stationary distribution. `finalize` was changed to write
  the report and `table2.csv` without a chain. That was a code fix after the
  pipeline stopped at this point; the run results are unchanged.
- **Fit** (OOS):
  - level −1.0%, var Δ1 −0.4%, var Δ5 +2.2%;
  - kurt Δ1 −4.6%, kurt Δ5 0.0%;
  - fractions +6.2% / −1.9% / −0.3%.
  - f = 0.0069, as good as v3. The model has no trouble fitting the
    end-of-window targets; it just wants no mean reversion.
- **Why β₂ → 0.**
  - The v3_end targets combine a lower level (0.631) with relatively high
    kurtosis and a var Δ5 that is low relative to the level.
  - In a lifecycle panel from z = 0, β₂ is identified only by the *curvature*
    of variance by age over 31 years. At β₂ ≈ 0 variance grows linearly:
    v3_end gives var at 25/35/45/55 = 0.11/0.46/0.80/1.14, against v3's
    0.13/0.59/0.87/1.05. So v3_end has lower dispersion at young ages and
    higher at 55.
  - Two of the targets that drive this (var Δ5, kurt Δ5) are **4-year
    extrapolations** (§26.4). The corner solution therefore rests partly on
    an extrapolated linear trend through a series with a level shift around
    2008–09.
- **Bound issue, decision.**
  - The bound is not lowered. Below 10⁻⁴ the half-life is already longer
    than the working life, so within the 31-year panel β₂ = 10⁻⁴ and β₂ = 0
    are observationally the same. Lowering the bound would only move the
    point along a flat direction.
  - v3_end stays a **robustness check only**, as the brief says. It is not
    usable as a HANK input as it stands: the stationary step-2 model needs an
    ergodic distribution, which this process does not have.
  - What it shows: **the end-of-window data push toward more, not less,
    persistence**. If anything, that strengthens the "high long-run
    dispersion" conclusion.

**v2 vs v3 vs v3_end**

| | v2 | **v3 (baseline)** | v3_end |
|---|---|---|---|
| λ₁ / β₁ / σ₁ | 0.162 / 2.83 / 2.00 | 0.162 / 2.80 / 2.00 | 0.158 / 2.82 / 1.95 |
| λ₂ / β₂ / σ₂ | 0.0099 / 0.0064 / 1.21 | 0.0100 / 0.0057 / 1.19 | 0.0068 / **0.0001 (bound)** / 1.13 |
| OOS f | 0.0073 | 0.0076 | 0.0069 |
| var₂ (additive formula) | 1.13 | 1.26 | none (random walk) |
| ergodic sd z₂ | 1.06 | 1.12 | none |
| half-life z₂ (yr) | 27 | 31 | ≈ 1700 |
| half-life z₁ (yr) | 0.061 | 0.062 | 0.061 |
| var(log y) at age 55, model | n/a (36-yr panel) | 1.05 | 1.14 |
| worst-fit moment | frac<10 +6.6% | frac<10 +6.4% | frac<10 +6.2% |

---

## 27. Asymptotic SMM standard errors for Poland v3 (C, 2026-11-08)

Script: `scripts/pol_v3_smm_se.py`. Output: `output/pol_v3/smm_se/{smm_se.csv,
smm_se.json}`. The estimate is the canonical v3 (§26.5).

    V = (G'WG)⁻¹ G'W (Ω_data + Ω_sim) W G (G'WG)⁻¹

**Ingredients.**
- **W = diag(1/m̂²).** These are exactly the weights of our objective. W is
  not the efficient weighting matrix, hence the full sandwich.
- **G = ∂m/∂θ** (8 × 6), by central differences in log θ on a 200k-worker
  panel with common random numbers (seed 1234, ages 25–55).
  - Two step sizes, h = 0.02 and 0.05 in logs.
  - Median elementwise relative difference between them: 0.45%. The largest
    differences are in single near-zero elements: σ₁ column 34%, β₂ column 23%.
  - **Resulting SEs change by ≤ 2% for β₂/var₂ and ≤ 10% for σ₂ and λ₂.**
    So G is stable enough. The table below uses h = 0.02.
- **Ω_data**: covariance of the annual statistics divided by the number of
  years, because each target is a mean over years. The series have different
  lengths (level 21, Δ1 20, Δ5 16 years), so element (j, k) is
  s_jk · n_jk / (T_j T_k), with s_jk taken over the n_jk common years. Two
  versions:
  - **detrended**: residuals from the §26.3 linear trends (divisor n − 2);
  - **raw**: deviations from the window mean (divisor n − 1).
- **Ω_sim**: covariance of the model moments across 60 simulation seeds, on
  the estimation-size panel (50k workers). One panel with CRN is used in
  estimation, so the simulation error enters once, with no 1/S factor.

SD of each moment, data (detrended / raw) vs simulation:

| | var log | var Δ1 | var Δ5 | kurt Δ1 | kurt Δ5 | frac<10 | frac<20 | frac<50 |
|---|---|---|---|---|---|---|---|---|
| data, detrended | 0.0070 | 0.0030 | 0.0040 | 0.237 | 0.106 | 0.0083 | 0.0051 | 0.0020 |
| data, raw | 0.0114 | 0.0059 | 0.0128 | 0.264 | 0.206 | 0.0087 | 0.0061 | 0.0027 |
| simulation (50k) | 0.0052 | 0.0010 | 0.0022 | 0.129 | 0.067 | 0.0005 | 0.0004 | 0.0004 |

The data SDs match the bootstrap target SDs in §26.6 to within 10–30%. That is
the expected agreement between an analytical s/√T and a resampling SD.
Simulation noise matters only for the level and the kurtoses, where it is
50–75% of the detrended data SD.

**Standard errors** (relative SE in %; "+sim" adds Ω_sim):

| | detrended + sim | detrended | raw + sim | raw | sim only | boot (ii) detrended | boot (i) years |
|---|---|---|---|---|---|---|---|
| λ₁ | 5.2 | 4.8 | 5.0 | 4.7 | 1.9 | 1.6 | 2.4 |
| β₁ | 6.0 | 5.6 | 5.9 | 5.4 | 2.3 | 2.1 | 2.4 |
| σ₁ | 0.5 | 0.4 | 0.7 | 0.6 | 0.3 | 0.4 | 0.5 |
| λ₂ | 3.2 | 2.8 | 4.1 | 3.8 | 1.5 | 2.5 | 4.2 |
| β₂ | **6.7** | 6.3 | **15.4** | 15.2 | 2.3 | **5.6** | **11.1** |
| σ₂ | 1.1 | 0.8 | 1.8 | 1.6 | 0.8 | 0.9 | 1.6 |
| var₂ (abs. SE) | **0.059** | 0.053 | **0.126** | 0.123 | 0.026 | **0.050** | **0.097** |
| erg. sd z₂ | 0.026 | 0.024 | 0.056 | 0.055 | 0.011 | 0.022 | 0.042 |
| half-life z₂ (yr) | 2.1 | 1.9 | 4.7 | 4.6 | 0.7 | 1.9 | 3.7 |

Derived quantities use the delta method with numerical gradients.

**Comparison with the bootstraps.**
- **The persistent component agrees.**
  - Asymptotic detrended vs bootstrap (ii): var₂ 0.053 vs 0.050, β₂ 6.3% vs
    5.6%, λ₂ 2.8% vs 2.5%, σ₂ 0.8% vs 0.9%.
  - Raw vs bootstrap (i): var₂ 0.123 vs 0.097, β₂ 15% vs 11%, σ₂ 1.6% vs
    1.6%, λ₂ 3.8% vs 4.2%.
  - So the conclusions in §26.6 do not depend on the method. The raw
    asymptotic SE is about 25% above the years bootstrap for β₂/var₂, which is
    plausible for a nonlinear map evaluated at T = 16–21.
  - With Ω_sim added (an honest total), var₂ = 1.26 ± 0.06 (detrended) or
    ± 0.13 (raw). That gives 95% intervals of about [1.14, 1.38] and
    [1.01, 1.51]. Both exclude the KMV reference (0.65, persistent
    component of the published chain, §29.9) by a wide margin. They also
    exclude the superseded 0.91.
- **The transitory component does not agree.** λ₁ and β₁ have asymptotic SEs
  of 5–6%, against 1.6–2.4% from the bootstrap. Two candidate reasons; the
  data here cannot separate them:
  - **The bootstrap understates.** The reduced budget (DE 10×5, NM 50, warm
    start) moves λ₁/β₁ only a little from the starting point. The validation
    in §26.6 shows reduced-vs-reference deviations of −2% on average for λ₁
    and β₁, the same order as the bootstrap SE itself. So the bootstrap does
    not resolve these parameters well.
  - **The asymptotic SE overstates.** The σ₁ column of G is the least stable
    (34%), and the information on (λ₁, β₁) comes almost entirely from the
    near-zero fractions and kurt Δ1. A linearization there is rough.
  - In either case the transitory parameters are well determined in economic
    terms (half-life 3 weeks ± a few days). They do not affect the step-2
    conclusions, which run through the persistent component.
- G'WG has condition number 3·10⁷. This reflects the scale of the parameters
  (β₂ ≈ 0.006 against β₁ ≈ 2.8), not near-singularity in logs. The relative SEs
  (`relse_*` in the CSV) are scale-free.

**Caveat (required by the brief).** Ω_data is **the variability of the annual
GRID statistics across years**, not sampling error from the microdata behind
each annual statistic. GRID publishes no microdata standard errors, and with
roughly 10⁵–10⁶ workers per year that error is probably smaller than the
year-to-year variation. These SEs therefore answer "how much would the estimate
move if the 2000–2020 years had been different draws from the same
(detrended or stationary) process", not "how precisely is a given year
measured". Two further assumptions:
- years are independent. The residuals carry the 2008–09 level shift, so they
  are not, and both the asymptotic and the bootstrap SEs are probably
  somewhat too small;
- the model is correctly specified. The J-type misfit in kurt Δ1 and
  frac < 10% (§26.5) is not reflected in the SEs.

---

## 28. Poland v3 vs KMV: figures, process table, grid size (D, E, F; 2026-11-08)

The Poland side uses the canonical v3 estimate (§26.5) and the KMV side the
Table 3 parameters. Continuous-process statistics are simulated with the same
code and design as the estimation: additive jumps, a lifecycle panel from
z = 0, ages 25–55. "KMV published chain" is
`annual-review/inputs/kmv_process` (3 × 11, the chain step 2 uses). Under §29.5
that chain was built for KMV's *reset* process, so it is not our process at KMV
parameters.

### 28.1 Figures (D)

Scripts and files:
- script `scripts/pol_v3_figures.py`;
- output `output/pol_v3/figures/`;
- numbers in `figure_stats.json`.

| file | content | reading |
|---|---|---|
| `d1_distribution.png` | 1-yr changes: CDF with the GRID percentiles (mean 2000–19), log density, one-sided log survival P(Δ > x) and P(Δ < −x) | see below |
| `d5_distribution.png` | the same for 5-yr changes. GRID gives no Δ5 percentiles, so the data enter as a normal with the target variance. | Poland and KMV are close; both are far from normal |
| `qq_vs_data.png` | model vs GRID quantiles of Δ1 (p1 … p99.99) | see below |
| `var_by_age.png` | var(log y) by age, 25–60 (56–60 shaded: outside the v3 window) | Poland and KMV coincide up to age ≈ 38; then Poland keeps rising, to 1.05 at 55 vs KMV 0.92. Mean over 25–55: 0.691 vs 0.661 |
| `sample_paths.png` | quarterly z₁, z₂, z₁ + z₂ for 6 workers | Poland's z₁ is a sequence of very short spikes. z₂ moves in rare jumps with slow drift |
| `impulse_response.png` | E[z(t)] after a +1 sd jump, both components | see below |
| `ergodic_chain.png` | ergodic distribution of log e on the chains | Poland 5×15 has var 1.51; KMV at Table 3 on our 5×15, 1.16; KMV published 3×11, 0.85 |

The tail panel was first drawn as P(|Δ| > x) against GRID's one-sided
percentiles, which overstates the model tail by about 2×. It was redrawn
one-sided before this section was written. The readings below refer to the
corrected figure.

**1-year changes vs GRID** (quantiles of Δ1):

| p | 1 | 5 | 25 | 37.5 | 50 | 62.5 | 75 | 95 | 99 | 99.9 | 99.99 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| GRID PL | −2.00 | −0.71 | −0.04 | +0.011 | +0.048 | +0.089 | +0.150 | +0.73 | +1.61 | +2.51 | +3.15 |
| Poland v3 | −1.68 | −0.75 | −0.11 | −0.024 | 0.000 | +0.019 | +0.106 | +0.77 | +1.71 | +3.04 | +4.29 |
| KMV T3 | −1.49 | −0.69 | −0.08 | −0.008 | 0.000 | +0.003 | +0.063 | +0.72 | +1.65 | +2.94 | +4.18 |

- **The body fits well**, from p5 to p99, within about 0.1 log points.
- **The centre is shifted.** GRID's residualized changes have median +0.05,
  while the model is symmetric around 0 by construction. This is the source of
  the "too many changes < 10%" misfit: the model's mass at 0 sits where the
  data have little, because the data are centred at +5%. This is a property of
  the target, not a parameter problem. The model has no drift term and cannot
  fit a non-zero median.
- **Tails are asymmetric in the data and symmetric in the model.**
  - The data's lower tail at p1 is heavier: −2.00 vs −1.68. That is large
    falls, i.e. exits to near-zero earnings.
  - The extreme upper tail is much thinner in the data: p99.99 is 3.15 vs
    4.29, and P(Δ > 2.5) is about 0.001 in GRID vs about 0.003 in the model.
  - So the model's kurtosis comes partly from upper-tail mass the data do not
    have. KMV's process has the same feature.
  - This cannot be fixed within the KMV family (symmetric normal jumps). It is
    a known limitation, reported here and not acted on.

**Impulse responses.**
- **Transitory component.** Poland: σ₁ = 2.00, half-life 0.062 yr (≈ 3 weeks).
  KMV: 1.74, 0.23 yr. Poland's "transitory" shocks are larger and almost
  immediately undone. Within a year they show up mostly as a mean shift of
  annual earnings, scaled by the share of the year left. KMV's last about a
  quarter.
- **Persistent component.**
  - Poland: σ₂ = 1.19, half-life 31 yr. KMV under the additive reading: 1.53,
    19 yr.
  - Poland's persistent jumps are smaller but much more durable. After about 18
    years the two responses cross.
  - Under KMV's own reset reading (§29.5) the expected response decays at
    β + λ, a half-life of ≈ 11 yr (dashed line). On that reading Poland's
    persistence is even further above KMV's.

### 28.2 Process table: Poland v3 vs KMV (E)

Script `scripts/pol_v3_vs_kmv.py table`. Output in
`output/pol_v3/vs_kmv/{vs_kmv.md, vs_kmv.json}`. Rates are quarterly.

**Parameters and economic quantities**

| | Poland v3 | KMV Table 3 |
|---|---|---|
| λ₁ / β₁ / σ₁ | 0.162 / 2.80 / 2.00 | 0.080 / 0.761 / 1.74 |
| λ₂ / β₂ / σ₂ | 0.0100 / 0.00567 / 1.19 | 0.007 / 0.009 / 1.53 |
| years between transitory shocks | 1.54 | 3.13 |
| years between persistent shocks | 24.9 | 35.7 |
| half-life z₁ / z₂ (yr) | 0.062 / 30.6 | 0.23 / 19.3 |
| ergodic sd z₁, additive formula λσ²/(2β) | 0.340 | 0.399 |
| ergodic sd z₂, additive formula | **1.122** | **0.954** |
| ergodic sd z₂, reset formula λσ²/(λ+2β) | 0.817 | 0.810 |
| ergodic var(z₁+z₂), additive / reset | 1.376 / 0.781 | 1.069 / 0.807 |

The reset-formula row applies KMV's specification to *both* parameter sets.
Poland's parameters were estimated under additive jumps, so "Poland, reset
formula" is **not** a valid Poland number. It only shows that the gap between
the two rows comes almost entirely from the formula.

**KMV reference (§29.9): the published chain** (3 × 11, bit-identical to the
replication package's discretization output):

| | var(log e) total | z₁ | z₂ |
|---|---|---|---|
| ergodic | **0.846** | 0.199 | **0.647** |

Against that, the Poland v3 chain (additive, 5 × 15) is 1.51 / 0.13 / 1.39.
The specification-consistent Poland number, re-estimated under resets, is in
§29.8.

**Continuous-process moments** (simulated, 100k workers, ages 25–55; in brackets
25–60)

| moment | Poland v3 | KMV Table 3 | PL target v3 |
|---|---|---|---|
| var log earns | 0.691 (0.746) | 0.650 (0.689) | 0.691 |
| var Δ1 | 0.262 | 0.227 | 0.260 |
| var Δ5 | 0.476 | 0.501 | 0.477 |
| kurt Δ1 | 14.41 | 14.61 | 15.24 |
| kurt Δ5 | 9.97 | 10.83 | 9.88 |
| frac < 10% | 0.479 | 0.548 | 0.450 |
| frac < 20% | 0.656 | 0.659 | 0.669 |
| frac < 50% | 0.836 | 0.840 | 0.843 |

- At ages 25–55 KMV's process gives a level of 0.65, below its own 0.70
  target. That target was matched on 25–60 under our design (0.69).
- Compared with KMV, Poland has more 1-year variance, less 5-year variance,
  and fewer small changes.

**Chains** (log e; annual statistics from the quarterly generator;
simulated chain moments on ages 25–55, 100k workers, starting from the state
nearest z = 0)

| | Poland v3 chain (5×15) | KMV T3, our 5×15 | KMV published (3×11) |
|---|---|---|---|
| states | 75 | 75 | 33 |
| ergodic var(log e) | **1.513** | 1.161 | 0.846 |
| P(stay in same state), annual | 0.849 | 0.803 | 0.846 |
| P(stay), quarterly | 0.901 | 0.902 | 0.925 |
| expected years in top decile | 3.12 | 3.05 | 1.83 |
| expected years in bottom decile | 3.12 | 3.05 | 1.83 |
| annual autocorrelation of log e | 0.894 | 0.827 | 0.726 |
| ergodic mass on the extreme states | 2e-5 | 1e-5 | 5e-4 |
| chain var log earns | 0.692 | 0.647 | 0.520 |
| chain var Δ1 | 0.175 | 0.218 | 0.240 |
| chain var Δ5 | 0.406 | 0.488 | 0.507 |
| chain kurt Δ1 | 17.0 | 14.1 | 11.5 |
| chain kurt Δ5 | 10.6 | 10.4 | 9.6 |
| chain frac < 10 / 20 / 50% | 0.66 / 0.68 / 0.83 | 0.65 / 0.73 / 0.87 | 0.72 / 0.73 / 0.82 |

The deciles are symmetric because both chains are symmetric in log e.

Readings:
- **Poland's chain is more persistent and more dispersed than either KMV
  chain.** Its autocorrelation is 0.89, against 0.83 and 0.73. Its ergodic
  variance is 1.51, against 1.16 and 0.85.
- **The two KMV chains differ a lot from each other.** Our discretization of
  KMV's parameters gives an ergodic variance of 1.16, while KMV's published
  chain gives 0.85. That gap is the additive vs reset specification (§29.5).
- For the HANK block, **the right benchmark for "Poland vs KMV" in step 2 is
  KMV T3 on our 5×15 chain**, which uses the same method and the same
  specification. Using the published 3×11 chain would add the specification
  difference to the country difference.
- Our 5×15 chain has very little mass on the extreme states, because the grid
  is wide. KMV's 3×11 has 25× more.

### 28.3 Grid size (F)

Script: `scripts/pol_v3_vs_kmv.py grid`, run as the **last** pipeline step,
with no other computation on the machine. Output: `output/pol_v3/vs_kmv/grid_size.*`.

Setup:
- the Poland v3 chain on four grids, each with `match_var_log` at ages 25–55;
- chain moments from 100k simulated workers, ages 25–55;
- the time is one steady-state solve of the `hank_step2` household block
  (`hh_ha.steady_state` at the baseline betas, median of 5 runs after a
  warm-up);
- relative gaps are chain vs continuous process.

| grid | states | ss solve (s) | ergodic var(log e) | var log | var Δ1 | var Δ5 | kurt Δ1 | kurt Δ5 | frac<10 | frac<20 | frac<50 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| KMV published 3×11 | 33 | 0.38 | 0.85 | | | | | | | | |
| 3×11 | 33 | 0.38 | 1.516 | +0.1% | −33.8% | −14.2% | **−22.1%** | −6.5% | +56.9% | +16.5% | +3.7% |
| 3×15 | 45 | 0.66 | 1.517 | +0.2% | −34.0% | −15.1% | −21.2% | −4.1% | +55.8% | +16.3% | +4.1% |
| 5×11 | 55 | 0.84 | 1.512 | +0.1% | −32.8% | −13.9% | **+16.9%** | +3.6% | +38.7% | +3.2% | −1.0% |
| **5×15 (default)** | 75 | 1.14 | 1.513 | +0.1% | −33.0% | −14.7% | +18.1% | +6.0% | +37.7% | +3.0% | −0.6% |

Readings:
- **The number of points for the transitory component (3 vs 5) matters; the
  number for the persistent one (11 vs 15) does not.**
  - With 3 points, kurt Δ1 is 22% too low and the < 10% fraction 56% too
    high.
  - With 5 points, kurt Δ1 is 17–18% too high and the < 10% fraction 38% too
    high.
  - 11 → 15 persistent points changes no moment by more than 2 pp.
- **var Δ1 is about −33% on every grid.** This is the §8d bias, a property of
  discretizing a fast-reverting component at quarterly frequency, not of grid
  size. Its fix stays on the open list (§8d).
- **Time is not a constraint.** Even 75 states solve in about 1.1 s.
  - An earlier smoke test under full CPU load (v2 parameters, while DE was
    running) gave 14.5 s and 122 s for the 5-point grids, and 0.57 s for
    the KMV 3×11 reference.
  - The 3×11 reference slowed only 1.5×, the 5-point grids about 100×, so
    contention alone does not explain the difference. The cause was not
    investigated: the parameters also differed (v2 vs v3). The clean final
    measurement above is the one to use, but if step 2 ever shows slow
    solves on 5-point chains, this is where to look.
  - The "chain size for step 2" item in Open decisions is therefore settled:
    **keep 5×15** (§14). If speed ever matters, 5×11 gives the same moments
    at 0.84 s.

---

## 29. How KMV (and GKOS) actually did it — and where we differ (G, 2026-11-08)

**Sources obtained.** Saved locally to `inputs/kmv_method/` and `inputs/gkos/`
(not committed: third-party PDFs and code):
- KMV, AER 108(3) 2018, published version (`HANK.pdf`), and the online
  appendix (`app.pdf`, benjaminmoll.com);
- the 2016 CEP working paper (`CEP-WP-2016-2-…pdf`);
- GKOS NBER w20913 (2015, the version KMV cite as "Guvenen et al. 2015"),
  GKOS 2019, and the moment workbooks `GKOS_2015_moments.xls` /
  `GKOS_2021_moments.xlsx` (github.com/serdarozkan/papers);
- `estimate.cu`, a CUDA port of "Kaplan, Moll and Violante's original Fortran
  codes" by Karibzhanov (Bank of Canada SDP 2020-13, `boc.pdf`;
  github.com/ikarib/HANK_EarningsEstimation).

**Could not be fetched; please download manually:**
1. **KMV replication package**, openICPSR project 231376
   (https://www.openicpsr.org/openicpsr/project/231376/version/V1/view, DOI
   10.3886/E231376V1, ~178 MB, needs a login). Unzip to
   `inputs/kmv_replication/`. This is the only primary source for points 1 and
   5 below; everything about the code is otherwise from the port.
2. **GKOS Econometrica 2021 + Supplemental Material** (ECTA14603):
   https://www.econometricsociety.org/publications/econometrica/2021/09/01/what-do-data-millions-us-workers-reveal-about-lifecycle,
   into `inputs/gkos/`. Low priority: KMV used the 2015 version.

### 29.1 The SMM objective

The paper and appendix do not define it. AER p. 719 says only "Simulated
Method of Moments using Social Security Administration (SSA) data". Appendix
D.1 covers only the discretization (3 × 11 points). The code port
(`estimate.cu`, lines 254–263 and 376):
- residuals are **relative**: `model/target − 1` for all 8 moments, as in ours
  (§4);
- **the two kurtosis residuals are multiplied by √0.5**, so in the squared
  objective the kurtoses get **weight 0.5**, and every other moment weight 1.
  (Our subagent's summary missed this; read directly from the file.) We use
  weight 1 for all 8 moments;
- an identity matrix otherwise. The optimizer is DFLS/NEWUOA (derivative-free
  least squares), not DE/TikTak.

### 29.2 Standard errors

None are reported. Table 4 (AER p. 720, "Rates expressed as quarterly
values") gives point estimates only; the appendix and the 2016 WP have none
either. **So our bootstrap (§19/§24.11/§26) and the asymptotic SEs (§27) go
beyond KMV**, and there is no KMV benchmark to compare them with.

### 29.3 The 0.70 target for the level

- **The seven change moments** come from GKOS w20913, Figure 1 / Table II:
  **one base year (1995–96 changes), pooled over all ages 25–60** of the base
  sample (men, above Y_min ≈ one quarter of full-time work at half the
  minimum wage, no large self-employment income). Changes are taken net of
  age effects. The subagent recomputed them from the workbook: var Δ1 0.234,
  kurt 17.80, var Δ5 0.459, kurt 11.55, fractions 0.53/0.71/0.86. The port's
  comment line gives the same targets (5-year kurtosis 11.55, not the 11.6 in
  `KMV_US_DATA_TARGETS`).
- **The 0.70 is defined nowhere**, neither in GKOS tables nor in KMV.
  - The nearest published object is w20913's Figure A.2, a within-cohort
    variance by age (0.59 at 25, 1.05 at 60; mean 0.79 over 25–60, 0.76 over
    25–55).
  - Most likely it is a cross-sectional variance of log earnings for the 1995
    base sample, possibly net of age dummies, supplied by GKOS and
    unpublished. Unverified; only the replication package, or an e-mail to
    the authors, can settle it.
- Consequence for us: there is no "correct" KMV definition to copy. v3's
  choice (age-residualized, men 25–55, mean over years) is defensible on its
  own terms. It is not "what KMV did", because that is unknown.

### 29.4 How GKOS aggregate over years

- **The KMV targets use a single year, pooled across ages**, not an average
  over years.
- Elsewhere GKOS average annual statistics: in w20913 the age/RE-cell moments
  are "computed for each year and then averaged over all years" (Table I
  note). In the 2019/2021 version, "average of these moments … over years".
- We (GRID-based targets) **average annual statistics over 2000–2020**. For
  kurtosis, pooling observations over years would add a mixture component
  whenever variance moves over time. Checked on the Polish year panel with the
  zero-mean mixture approximation, using nobs weights:
  pooled kurt Δ1 ≈ 15.20 vs the mean of annual values 15.24, and pooled kurt
  Δ5 ≈ 9.84 vs 9.88. **Negligible**: over this window the variance falls
  while kurtosis rises, so the mixture effect cancels. No change needed.

### 29.5 The jump specification: KMV's published process uses RESET jumps (new finding, important)

§1 decided "additive jumps, not resets", citing "KMV eq. 30–31":
`dz = −βz dt + ε dN`. **Those equation numbers are from the 2016 CEP working
paper.** The **published AER version** (p. 718–719, eqs. 33–34) says instead
that conditional on a jump, "a new log-earnings state z_{j,it} is drawn from
a normal distribution" N(0, σ²), with `dz = −βz dt + dJ`. That is a **reset**:
on arrival the component is replaced, not incremented. The online appendix
(`app.pdf`, the y-process paragraph) uses the same "new log-earnings state is
drawn" wording. Two independent checks agree with the published text:

1. **KMV's own chain** (`annual-review/inputs/kmv_process/ymarkov_combined.txt`,
   the generator used in step 2). With additive jumps the rate of jumping
   *into* state j depends on the origin state (a normal centred on the
   origin). With resets it is λ·w_j for every origin. In the 33-state
   generator there are only **13 distinct off-diagonal rates**, and at most 3
   distinct rates flow into any state (persistent-jump landings 0.00036–0.0008,
   identical across origins, plus drift entries). That is the reset signature.
   Script: `scripts/kmv_port_check.py` → `output/pol_v3/kmv_method/port_check.json`.
2. **The code port**: `z = jump ? σ·ε : (1 − δ/4)·z` (reset), 16 steps/year,
   start from N(0, σ²/(1 + 2δ/λ)) (the reset process's stationary variance),
   25-year burn-in, then moments in **one year of an ergodic cross-section**,
   not on a lifecycle panel. Its "5-year" change is `log(y_4/y_0)`, i.e. a
   4-year gap (probably a bug in the port, or in the original).

**Does the reset + ergodic design reproduce KMV's Table 2?** Simulated at Table 3
parameters (400k workers):

| moment | KMV Table 2, model | ours: additive + lifecycle 36 yrs (§1–2) | port: reset + ergodic, 4-yr "Δ5" | reset + ergodic, true 5-yr Δ5 |
|---|---|---|---|---|
| var log earns | 0.70 | 0.693 | 0.750 | 0.750 |
| var Δ1 | 0.23 | 0.229 | 0.234 | 0.234 |
| var Δ5 | 0.46 | 0.504 | 0.484 | 0.545 |
| kurt Δ1 | 16.5 | 14.55 | **16.18** | 16.18 |
| kurt Δ5 | 12.1 | 10.84 | 12.95 | **11.90** |
| frac < 10% | 0.56 | 0.547 | 0.552 | 0.552 |
| frac < 20% | 0.67 | 0.659 | 0.663 | 0.663 |
| frac < 50% | 0.85 | 0.839 | 0.845 | 0.845 |

- **Both designs reproduce Table 2 about equally well overall.** Ours is
  better on the level and var Δ5. Reset + ergodic is clearly better on both
  kurtoses, which is exactly the "residual kurtosis gap" the `simulate.py`
  docstring attributes to GKOS sample details.
- So **§1's verification ("only the additive version reproduces"; reset gives
  var Δ1 ≈ 0.40) does not hold for the reset + ergodic design**. It was
  presumably run as reset + lifecycle, or with a different reset
  implementation. That test is not in the repo and cannot be re-checked.
- Our design is internally consistent and reproduces Table 2. **But it is the
  2016 WP's process, not the published one.** The published process, KMV's own
  chain and the port all use resets.

**A third check: the ergodic variance of KMV's chain.** The published 33-state
chain has ergodic var(log e) = **0.846** (`hank_step2/results/baseline.json`).
At Table 3 parameters the reset process implies λσ²/(λ + 2β) summed over the
two components = 0.151 + 0.655 = **0.81**; the additive process implies
**1.07**. STEP2_LOG §S5 read the 0.85 as "a specific 33-state discretization,
not KMV's target" and took 1.07 as KMV's value. Under the published (reset)
specification 0.85 is simply the chain matching its own process to within 5%.
That 1.07 anchor (STEP2_LOG, §25 table, the GER/USA pin) is the additive formula
applied to reset estimates.

**Why it matters for the thesis conclusion.**
- Under resets, the ergodic variance of a component is λσ²/(λ + 2β), not
  λσ²/(2β). The "KMV anchor" used throughout (§8b, `KMV_ERGODIC_VAR2 = 0.9103`
  for the persistent component) applies the **additive** formula to
  parameters KMV estimated under **resets**. Under KMV's own specification
  their persistent component has ergodic variance 0.007·1.53²/(0.007 + 0.018)
  = **0.655** (sd 0.81), not 0.91 (sd 0.95).
- Comparisons of Poland's var₂ with "KMV's value" (§23, §24.5, §24.12)
  therefore mix two specifications. The like-for-like comparison is either
  Poland-additive vs KMV-re-estimated-additive (§12's GRID-USA is the nearest
  thing), or Poland-reset vs KMV-reset.
- β has a different meaning under the two specifications. It is the drift
  between jumps in both, but under resets a jump also wipes out past
  persistence. Half-lives and the "years between shocks" are not comparable
  across the two.

**Decision (2026-11-08).** The B–F runs started today use the **existing
design** (additive jumps, lifecycle panel), so v3 stays comparable with v1/v2
and with FRA/GER/USA, which the brief says not to recompute. All "vs KMV"
comparisons in §26–§28 are labelled with the formula used (additive-formula
anchor 0.91 *and* KMV-spec 0.655). **Recommended follow-up, not run:** add a
`reset_jumps` option to the simulator (default unchanged). Re-estimate Poland
v3 under resets, with a lifecycle panel and with an ergodic cross-section,
and decide with the supervisor which specification is the baseline. This
needs the openICPSR package first (point 1 above), to confirm the original
Fortran matches the port.

### 29.6 Summary: our method vs KMV

| | KMV (published / port) | ours | needs a change? |
|---|---|---|---|
| objective | squared relative deviations, identity, **kurtosis weight 0.5** | squared relative deviations, all weights 1 | **minor**: a weight-0.5 variant is a cheap robustness check |
| optimizer | DFLS / NEWUOA | DE + TikTak multistart, OOS selection | no (ours is more thorough) |
| SEs | none | bootstrap + asymptotic (§27) | no |
| jump specification | **reset** (AER eq. 33–34, KMV chain, port) | additive (2016 WP eq. 31) | **yes, the most important one**: see 29.5 |
| simulation sample | ergodic cross-section, one year | lifecycle panel from z = 0, 36 (v3: 31) years | tied to the jump choice; decide together |
| time step | 16 / year, Euler decay | 24 / year, exact exponential decay | no |
| 5-year change | 4-year gap in the port (likely a bug) | true 5-year gap | no |
| data aggregation | one year (1995), pooled over ages | mean of annual statistics, 2000–2020 | no (pooling vs averaging is ≈ 0 for PL kurtosis, 29.4) |
| level target | 0.70, definition unknown | `std_resid²`, men 25–55 (v3) | no (KMV's is undocumented) |

### 29.7 The KMV replication package: what the original code does (2026-11-09)

**Source.** openICPSR 231376, unpacked locally at `data/KMV_replication/`
and not in git (`.gitignore`, commit 2666def). Below, `R/` means
`data/KMV_replication/231376-V1 (2)/20160042_data/Replication/`. File and line
numbers refer to that tree. No code from it is copied into this repository.
The two parts used are:
- `R/EarningsProcess/EarningsEstimation/` (SMM estimation, Fortran);
- `R/EarningsProcess/EarningsDiscretization/` (the 3 × 11 chain).

**The process in the code is a reset with a general ρ, and ρ is set to 0.**
- **On a jump** (`EarningsEstimation/Simulate.f90`, l. 43–47 and 55–59) the
  new value of the component is ρ · (old value) + σ · ε.
- **Between jumps** (l. 49 and 61, `AdditiveDrift = 0`) the value is
  multiplied by (1 − dt · β), i.e. Euler decay.
- **ρ is never estimated.** `SetParameters.f90` l. 30–31 set `EstimateRho1 =
  EstimateRho2 = 0`, and `MakeGuess.f90` l. 10 and 17 fix both ρ guesses at
  0.0. The paper's output `earnings_estimation_paper/parameters.txt` reports
  `rho1 = rho2 = 0`.
- So the estimated process is **z → ε on arrival: a reset**. An additive
  process would be ρ = 1.
- The stationary variance used for the starting draw (`Simulate.f90`
  l. 14/21) is λσ² / (2β + λ(1 − ρ²)). This is the reset formula at ρ = 0 and
  the additive one at ρ = 1.
- **The chain is built the same way.**
  - `EarningsDiscretization/TransitionMatrix.f90` l. 79–83 and 104–108: the
    landing probabilities are normal CDF differences centred at ρ ·
    (origin), with `rho1 = rho2 = 0.0` set in
    `EarningsDiscretization/SetParameters.f90` l. 26 and 31. So the jumps
    land on N(0, σ²) from every origin.
  - The drift is a 2-point linear-interpolation scheme
    (`DriftPointsVersion = 2`, `EarningsDiscretization/SetParameters.f90`
    l. 13).
  - The grid width and power spacing of each component are chosen by the
    same derivative-free least-squares routine (DFLS/NEWUOA) to match the
    same 8 moments, simulated on the chain.
  - `earnings_discretization_paper/ygrid_combined.txt` and
    `ymarkov_combined.txt` are **bit-identical** to
    `annual-review/inputs/kmv_process/` (checked: max abs. difference 0).
    That is the chain inside KMV's model and inside our S1 baseline.

**Estimation details** (what §29.1–29.4 inferred from the CUDA port, now
confirmed against the original):

| item | original code | where |
|---|---|---|
| objective | DFLS on residuals √w_m · (model/target − 1); w = 0.5 for kurt Δ1 and kurt Δ5, 1 otherwise | `EarningsEstimation/dfovec.f90` l. 126–185 (weights set at l. 145, 162) |
| optimizer | NEWUOA-based DFLS, 3 restarts (`ndfls = 3`), rhobeg 5, rhoend 1e-4 | `SetParameters.f90` l. 61–63; `Estimation.f90` l. 107–110 |
| parameter bounds | logistic transforms: λ ≤ 2, σ ≤ 2, β ≤ 1 (quarterly) | `SetParameters.f90` l. 65–69; `dfovec.f90` l. 25, 46, 82, 89, 96, 103 |
| simulation | dt = 0.25 quarter (16 steps/yr); start from the stationary normal; 100-quarter burn-in; 20 quarters kept = 5 annual observations; **one cross-section** | `Parameters.f90` l. 5–8; `Simulate.f90` l. 31–34, 73–78 |
| simulated workers | `nsim = 5000` as shipped (a comment says 50000; the discretization uses 50000) | `EarningsEstimation/Parameters.f90` l. 5; `EarningsDiscretization/Parameters.f90` l. 7 |
| "5-year" change | year 5 − year 1 = **a 4-year gap** (the bug suspected in the port is in the original) | `ComputeMoments.f90` l. 55–58 |
| targets | hard-coded: 0.7, 0.23, 17.8, 0.46, 11.55, 0.54, 0.71, 0.86; no data code in the package | `SetParameters.f90` l. 49–59 |
| standard errors | **none**: no Jacobian, Hessian or covariance anywhere in the estimation code | whole `EarningsEstimation/` |
| unrounded estimates | λ₁ 0.07987, β₁ 0.7616, σ₁ 1.7352, λ₂ 0.006561, β₂ 0.009392, σ₂ 1.5259 | `earnings_estimation_paper/parameters.txt` |

For §29.3: the package does not define the 0.70 either. It is a hard-coded
number with no data code behind it. That conclusion stands.

**How the paper describes the process.**
- The NBER working paper **w21897, section 4.3** (as noted by the user; the
  2016 CEP version we have, eq. 31, has the same form) writes the process as
  **dz = −βz dt + ε dN, i.e. additively**.
- The AER text (eq. 33–34, "a new log-earnings state is drawn") and the code
  are **reset**.
- So **in KMV the working paper's equation and the implementation differ.**
  The published paper and the code agree with each other.
- **Our baseline stays with the specification written in the working paper
  (additive).** The reset version is a robustness check (§29.8), not a new
  baseline.

**The contradiction with §1, resolved.** §1 (project start) found that
"reset at Table 3 parameters gives var Δ1 ≈ 0.40 against 0.23" and chose
additive jumps on that basis. The §1 test is not in the repository.

All designs were re-run at Table 3 parameters (`scripts/reset_vs_sec1.py`,
200k workers, `output/pol_v3/kmv_method/reset_vs_sec1.json`):

| design | var log | var Δ1 | var Δ5 | kurt Δ1 | kurt Δ5 | <10% | <20% | <50% |
|---|---|---|---|---|---|---|---|---|
| **KMV Table 2, model** | 0.70 | 0.23 | 0.46 | 16.5 | 12.1 | 0.56 | 0.67 | 0.85 |
| KMV code output (`moments.txt`, nsim 5000) | 0.696 | 0.229 | 0.464 | 16.53 | 12.06 | 0.557 | 0.673 | 0.851 |
| **D: replica of KMV's code**, unrounded params | **0.702** | **0.227** | **0.463** | 15.67 | 13.02 | 0.556 | 0.666 | 0.847 |
| E: replica, Table 3 rounded params | 0.754 | 0.233 | 0.486 | 15.72 | 12.90 | 0.554 | 0.663 | 0.844 |
| **A: ours** (additive, lifecycle 36 yrs from 0, true Δ5) | 0.691 | 0.228 | 0.504 | 14.53 | 10.79 | 0.547 | 0.659 | 0.840 |
| B: reset in our design (lifecycle from 0) | 0.573 | 0.227 | 0.520 | 15.70 | 11.75 | 0.551 | 0.661 | 0.842 |
| C: reset, ergodic, true 5-yr Δ5 | 0.751 | 0.232 | 0.555 | 16.01 | 11.87 | 0.545 | 0.657 | 0.841 |
| F: additive, ergodic, true 5-yr Δ5 | 1.017 | 0.230 | 0.520 | 14.47 | 10.55 | 0.540 | 0.655 | 0.838 |
| reset with all rates ×4 (units slip) | 0.684 | **0.421** | 1.061 | 10.45 | 5.72 | 0.276 | – | – |

The replica (D) follows the description above: reset, Euler decay, 16 steps
per year, stationary start, 25-yr burn-in, 4-yr gap, one cross-section. It
was written independently from the description, not copied. Last row: 60k
workers, run separately; < 20% and < 50% were not recorded.

**Comparing the two implementations of the reset.**
- **What is drawn.**
  - KMV: per step, a Bernoulli jump indicator with probability dt · λ, and on
    a jump a new value σ · ε (ρ = 0); the old value is discarded.
  - Our `reset_jumps=True`: the same, with probability 1 − e^{−λ dt}.
- **Decay in the same step.**
  - KMV: decay and jump are exclusive. A jumping step sets z = σε, a
    non-jumping step decays.
  - Ours: the same, because the decayed value is computed but unused on a
    jump. With additive jumps ours decays *and* adds in the jumping step.
- **Time step.** KMV uses 16 steps per year with Euler (1 − dt β). We use 24
  per year with exact e^{−β dt}. For β₁ = 0.76 per quarter the Euler step
  decays 19% per step against 17% exactly. This is a small effect on
  transitory persistence.
- **Sample.**
  - KMV: one ergodic cross-section, 5 annual observations, Δ5 over a 4-year
    gap.
  - Ours: a lifecycle panel from z = 0, 36 (v3: 31) years, a true 5-year gap.
- **Objective.**
  - KMV: kurtosis weight 0.5, DFLS with 3 restarts, 5000 simulated workers
    and one seed.
  - Ours: weight 1 for every moment, DE/TikTak, 50k workers, CRN.

**What explains §1's 0.40.**
1. **No correct reset implementation gives 0.40.** At Table 3 parameters
   every reset variant gives var Δ1 = 0.227–0.233 (B, C, D, E), exactly like
   the additive one. var Δ1 is driven by the transitory component, whose
   one-year change is nearly the same under reset and additive jumps (it
   decays within a quarter either way).
2. §1's number **is** reproduced by a reset with **all rates four times too
   fast** (var Δ1 0.42). That is quarterly rates applied per quarter-step when
   the step was meant to be a year or a sub-step. It is the most likely
   explanation: a units slip in the original check. It cannot be confirmed
   because the §1 script was not kept.
3. **What really separates the specifications is the level and Δ5, through
   the sampling design, not the jump type.**
   - In our lifecycle-from-0 design a reset gives a level of 0.57, too low.
     Persistent jumps cannot accumulate because each new one replaces the
     last.
   - In KMV's ergodic design a reset gives 0.70, as it should, while the
     additive process gives 1.02, too high.
   - Each specification reproduces Table 2 only in its own design: additive +
     lifecycle (A), or reset + ergodic + 4-yr gap (D).
   - §1 and §2 were decided one after the other, and their combination (A)
     reproduces Table 2. But the reason §1 gave ("reset misses var Δ1") is
     wrong.
4. **The replica reproduces KMV's own output at the unrounded parameters.**
   - The level and variances match to 0.001–0.006.
   - Kurtoses are 15.7 vs 16.5 and 13.0 vs 12.1. With nsim = 5000 the SD of
     the sample kurtosis is large: over 40 replica seeds at nsim = 5000 the SD
     is 3.0 for kurt Δ1 and 0.96 for kurt Δ5, mean 16.5 and 13.0. KMV's
     reported kurtoses are therefore single noisy draws, and the weight 0.5
     partly offsets that noise.
   - Rounding to Table 3 (λ₂ 0.00656 → 0.007, β₂ 0.00939 → 0.009) raises the
     level from 0.70 to 0.75 (E vs D). Table 3's rounded numbers do not quite
     reproduce Table 2 under KMV's own design. This also explains the 0.750 in
     §29.5's table, which used the rounded parameters.

**Corrections to §29.1–29.6 and §1.**
- **§29.5.** "the published process uses resets" is now confirmed in the
  original code, not just the port. The KMV reset variance is **0.603**
  (unrounded parameters) instead of 0.655 (rounded). §29.9 replaces both with
  the published chain.
- **§29.6.** The "4-year gap (likely a bug)" is in the original.
- **§1.** Its verification sentence is wrong, as shown above. The choice of
  additive jumps stands, on the ground of the working-paper specification
  (decision above), not on that test.

### 29.8 Poland v3 re-estimated with reset jumps (robustness, not a baseline) (2026-11-09)

**Setup.** `scripts/pol_v3_estimate.py --version v3_reset`, output
`output/pol_v3_reset/`.
- Identical to v3 (§26.5) except `reset_jumps=True`: the same targets
  (`targets/pol_grid_v3.json`), ages 25–55, lifecycle panel from z = 0,
  equal weights, CRN, the 5×15 chain.
- Starts: 3 DE seeds (7, 101, 202) + TikTak, no bootstrap, as briefed.
- Seed 202 was interrupted once (a broken pipe at 11:49, the machine
  sleeping) and restarted. Finished starts are skipped on restart, so nothing
  else was affected. The crash log is `logs/de_crash1.log`.

**All four starts land on one point, with β₂ at its lower bound.**
- The three DE seeds agree to 6 digits (OOS f 0.022525–0.022526). TikTak is
  0.5% worse (0.02270) at the same point.
- Canonical is `de_seed7`. β₂ = 1.0e−4 per quarter is the lower bound
  (log position 0.00). The persistent component therefore has **no drift**:
  it is constant between jumps and redrawn from N(0, σ₂²) on a jump.

| parameter / derived | v3 additive (§26.5, DE 404) | **v3 reset** (DE 7) | spread across reset starts |
|---|---|---|---|
| λ₁ | 0.1623 | 0.1772 | 0.2% |
| β₁ | 2.801 | 3.082 | 0.6% |
| σ₁ | 1.998 | 1.944 | 0.05% |
| λ₂ | 0.01004 | 0.01065 | 0.4% |
| β₂ | 0.00567 | **0.0001 (bound)** | – |
| σ₂ | 1.193 | 1.118 | 0.5% |
| var₂, formula (additive λσ²/2β; reset λσ²/(λ+2β)) | 1.260 | 1.228 | 1.1% |
| var₁, formula | 0.116 | 0.106 | |
| **OOS f** | **0.00762** | **0.02253** (×3.0) | 0.0225–0.0227 |
| in-sample f | 0.00734 | 0.02302 | |

**Fit (OOS relative errors, model/data − 1).**

| moment | data | v3 additive | v3 reset |
|---|---|---|---|
| var log | 0.691 | −0.1% | **−6.3%** (0.648) |
| var Δ1 | 0.260 | +0.8% | −2.4% |
| var Δ5 | 0.477 | −0.2% | **+8.8%** (0.519) |
| kurt Δ1 | 15.24 | −5.4% | −7.1% |
| kurt Δ5 | 9.88 | +1.0% | +4.2% |
| < 10% | 0.450 | +6.4% | +5.4% |
| < 20% | 0.669 | −2.0% | −1.8% |
| < 50% | 0.843 | −0.9% | −0.9% |

- The whole loss is in the **level and Δ5**. The transitory block (Δ1 moments)
  is fitted about as well as under additive jumps.
- This is §29.7 point 3 again, now in estimation. In a lifecycle panel from
  z = 0 a reset cannot accumulate persistent variance: each jump replaces the
  last one. Lowering β₂ to 0 is the only way left to raise the level. It
  still falls 6% short, and it overshoots Δ5, because the remaining way to
  build dispersion is large rare jumps.
- So reset + our design is a **binding misspecification**: it fits 3× worse,
  with a parameter on its bound. This is not evidence on which jump type is
  "true". Under KMV's ergodic design the ranking flips (§29.7, D vs F). It is
  evidence that the additive baseline is the internally consistent choice
  for our lifecycle targets.

**Ergodic variance of the chain (5×15, `match_var_log`).**

| | v3 additive | v3 reset | KMV published chain (§29.9) |
|---|---|---|---|
| **total var(log e)** | **1.51** | **1.43** | **0.846** |
| persistent z₂ (chain) | 1.39 | 1.32 | 0.647 |
| transitory z₁ (chain) | 0.127 | 0.114 | 0.199 |
| z₂ chain / formula | +10% | +7% | +7% (reset formula, unrounded) |

- Under resets Poland's long-run dispersion falls only 5% (total) and 5%
  (z₂). It is still **1.7× KMV's chain** in total and **2.0×** for the
  persistent component.
- So the main conclusion of v2/v3 (Poland is far more dispersed in the long
  run than KMV) **does not depend on the jump specification**. This is the
  like-for-like statement §30 asked for: reset Poland vs reset KMV, chain vs
  chain.
- Discretized moments (`table2.csv`) have the same pattern as in v3: var Δ1
  0.170 vs 0.253 and < 10% 0.67 vs 0.48. The 5-point transitory grid cannot
  reproduce the very short transitory spells (§28.3); the level, 0.646, is
  matched by construction.

**Half-lives of the persistent component.**

| | v3 additive | v3 reset |
|---|---|---|
| drift only, ln2/β₂ | 30.6 yr | 1733 yr (β₂ at the bound: no drift) |
| expected value, ln2/(β₂ + λ₂) (reset); ln2/β₂ (additive, jumps are mean zero) | **30.6 yr** | **16.1 yr** |
| mean years between persistent jumps, 1/λ₂ | 24.9 yr | 23.5 yr |
| transitory, expected | 0.062 yr | 0.053 yr |

- Under resets the drift half-life is meaningless (no drift). The relevant
  persistence measure is the autocorrelation of z₂, exp(−(β₂ + λ₂)t). Its
  half-life is **16 years**: a persistent state lasts until the next jump,
  on average 23.5 years, and is then forgotten entirely.
- Under additive jumps a jump only adds to the past, so the half-life is the
  drift's, 31 years.
- The "persistence" of Poland therefore halves under resets, while the
  dispersion barely moves. Step 2 should be read the same way: risk is
  about as large, but less permanent.
- Compared with KMV: KMV's z₂ autocorrelation half-life is
  ln2/(β₂ + λ₂) = ln2/(0.00939 + 0.00656)/4 = 10.9 yr. Poland reset is 16 yr,
  so it is still more persistent.

**Not done.** There is no bootstrap (as briefed) and no ergodic-design
estimation of Poland (not in the brief; it would also change how the targets
map to the model, cf. the 4-yr gap in §29.7).

### 29.9 The KMV reference point: the published chain, not a formula (2026-11-09)

**Decision.**
- Every "vs KMV" comparison of long-run dispersion now uses the **ergodic
  distribution of the published KMV chain**, `annual-review/inputs/kmv_process`.
  It is bit-identical to `R/EarningsProcess/EarningsDiscretization/
  earnings_discretization_paper/` (§29.7).
- This is the chain that sits inside KMV's model and inside our S1 HANK
  baseline. It is therefore what "KMV's earnings risk" means for the step-2
  comparison, whatever formula one would apply to their parameters.

| object | value | status |
|---|---|---|
| **published chain, ergodic var(log e)** | **0.846** (sd 0.920) | **reference** |
| published chain, persistent component z₂ (11 states) | **0.647** | reference for var₂ comparisons |
| published chain, transitory component z₁ (3 states) | 0.199 | |
| reset formula λσ²/(λ+2β), unrounded parameters: z₂ / z₁ | 0.603 / 0.150 | the continuous process behind the chain; the chain is +7% / +33% |
| reset formula at rounded Table 3 parameters (z₂) | 0.655 | superseded (§29.5) |
| additive formula at Table 3 parameters (z₂), `KMV_ERGODIC_VAR2` | 0.910 | superseded: our formula on their parameters, never in their model |
| additive formula, total | 1.07 | superseded (STEP2_LOG, §25, the GER/USA pin) |

**Places changed in this stage:** §26.5, §26.6, §26.7, §27, §28.2 and §30.
Each now quotes 0.85 (total) or 0.65 (persistent) and calls 0.91 "superseded".

**Places not changed:** §8b, §12, §23 and §24. They are earlier stages and
are kept as written, with this section as the erratum. The GER and USA
`β₂` pins (`ErgodicVarConstraint`, var₂ = 0.9103) were calibrated to the
superseded additive anchor.
- The brief says not to recompute FRA/GER/USA, so the pins stay.
- **This is now an open item**: they were pinned to a number KMV's model
  never contained. Re-pinned to the chain's 0.647 they would be more
  mean-reverting.

**Comparing like with like.**
- Chain vs chain, total var(log e):
  - Poland v3 additive (5×15): 1.51;
  - Poland v3 reset (5×15, §29.8): 1.43;
  - KMV published (3×11): 0.85.
- Persistent component, chain vs chain: Poland v3 1.39 (reset 1.32) vs KMV 0.65.
- The chain/formula ratio is similar for both (Poland +10%, KMV +7% on z₂).
  So formula-based Poland intervals (§26.6, §27) can be read against 0.65.

---

## 30. Summary of the v3 stage (2026-11-09)

**Done** (one commit per step, all pushed):

| item | where | result in one line |
|---|---|---|
| A1 v3 targets | §26.1, `targets/pol_grid_v3.json` | level = mean `std_resid²` 2000–20 = 0.691 (v2: 0.722) |
| A2 ages 25–55 option | §26.2, `tests/test_age_window.py` | `n_age_cols=31`; default bit-identical (test) |
| A3 trends | §26.3, `output/pol_v3/trends/` | variances fall significantly; kurt Δ5 rises; level shift around 2008–09 |
| A4 v3_end targets | §26.4, `targets/pol_grid_v3_end.json` | trend at 2019; var Δ5 and kurt Δ5 extrapolated 4 yrs |
| B1 estimation | §26.5 | 5 DE + TikTak, one basin; canonical DE 404, var₂ 1.26, hl₂ 31 yr |
| B2 v3_end | §26.8 | β₂ at the lower bound in every start (random walk); robustness only |
| B3 profile / recovery / nested | §26.7 | minimum at 1.0×; KMV's published chain (z₂ 0.65) below the lowest profiled point (0.756, f ×2.9); the old 0.91 already on the steep side; recovery 4/4; nested specs as in v2 |
| B4 bootstrap | §26.6 | var₂ 95%: [1.06, 1.44] (years), [1.17, 1.35] (detrended) |
| C asymptotic SEs | §27 | match the bootstraps for the persistent component; larger for λ₁/β₁ (recovery supports the asymptotic ones) |
| D figures | §28.1, `output/pol_v3/figures/` | body of Δ1 fits well; data centred at +5%, asymmetric tails |
| E vs KMV table | §28.2, `output/pol_v3/vs_kmv/` | Poland more persistent and more dispersed than KMV on every measure |
| F grid size | §28.3 | 3 vs 5 transitory points matters, 11 vs 15 does not; keep 5×15 (1.1 s) |
| G KMV method | §29.1–29.7 | replication package: KMV's code uses reset jumps (ρ = 0); kurtosis weight 0.5; no SEs; 0.70 hard-coded; §1's 0.40 was a units slip |
| reset robustness | §29.8, `output/pol_v3_reset/` | β₂ at the bound, OOS f ×3 worse; chain var(log e) 1.43 vs 1.51 additive; hl₂ 16 vs 31 yr |
| KMV anchor | §29.9 | published chain: total 0.846, z₂ 0.647 (replaces 0.91 / 0.655) |

**Does v3 change the main v2 conclusion (high long-run dispersion in Poland,
above KMV)? No. It strengthens it, with one important qualification.**
- v3's var₂ is 1.26 (sd 1.12), up from v2's 1.13 (sd 1.06), even though the
  level target fell 4%. The age window explains it (§26.5).
- Both bootstrap intervals and both asymptotic intervals exclude the KMV
  reference. That reference is now the published chain (§29.9): z₂ 0.647,
  total 0.846. They already excluded the superseded 0.91. v3_end goes
  further, to no mean reversion at all.
- **The specification question (§29.5, §29.7–29.8) is settled for this
  conclusion.**
  - KMV's code uses reset jumps; our baseline stays additive, following the
    working paper.
  - Re-estimated with resets, Poland's chain has var(log e) 1.43 (additive
    1.51) and z₂ 1.32 (1.39).
  - Against KMV's published chain (0.846 / 0.647) Poland is 1.7× / 2.0× in
    either specification, chain vs chain.
  - What the specification does change is **persistence**: the half-life of
    z₂ is 31 yr additive vs 16 yr reset (KMV reset: 11 yr). The reset fit is
    3× worse, with β₂ on its bound, so the additive baseline is also the
    better-fitting one in our design.

**Left / recommended next** (not run):
1. ~~Reset-jump re-estimation of Poland v3~~ **Done (§29.7–29.8)**: option
   `reset_jumps` (default bit-identical, `tests/test_reset_jumps.py`),
   lifecycle design, 3 DE + TikTak. Not run: an ergodic-design estimation and
   a reset bootstrap. Supervisor: confirm the additive baseline (§29.7).
2. A kurtosis-weight-0.5 variant of the objective (§29.6), a cheap robustness
   check.
3. Step 2 on the v3 chain (`output/pol_v3/income_process_*`):
   - redo the §24.6 one-/two-asset trigger check (v2 failed it; v3 is *more*
     dispersed, so a pass is unlikely);
   - the profile chains `output/pol_v3/var2_profile/chain_0.80`–`0.90` are the
     candidates for the lower-var₂ route inside the statistical band.
4. The fraction-construction gap (§26.1) and the 2008–09 level shift
   (§26.3/§26.6) are documented and not corrected.

---

## 31. Poland estimated with "the KMV method" (2026-11-10)

**Goal.** An estimation variant that is fully comparable with KMV's own
procedure, next to the v3 baseline (§26.5), which is not changed. Everything
below is derived from the replication package (openICPSR 231376, local only,
never committed). Only Poland is computed. Wherever the model has ages, they
match the data (Poland: 25–55). New outputs go to `output/pol_v3_kmvmethod/`.

**Path conventions.** `R/` = `data/kmv_replication/231376-V1 (2)/20160042_data/Replication/`,
`E/` = `R/EarningsProcess/EarningsEstimation/`, `D/` = `R/EarningsProcess/EarningsDiscretization/`.
"l." = line numbers in those files. Our code is cited as `kmv_earnings/…` and
`scripts/pol_v3_estimate.py` at commit 713f29f.

### 31.1 Inventory of differences: KMV's code vs our v3

**A. Process**

| # | element | KMV (file:line) | our v3 | differs? |
|---|---|---|---|---|
| A1 | jump type | on a jump `y ← ρ·y + σ·ε` (`E/Simulate.f90` l. 45–47, 57–59); ρ not estimated (`E/SetParameters.f90` l. 30–31), fixed at 0 (`E/MakeGuess.f90` l. 10, 17); reported ρ₁ = ρ₂ = 0 (`E/earnings_estimation_paper/parameters.txt` l. 7–8). So **reset**, z → ε | **additive**, z → z·e^{−βdt} + ε (`kmv_earnings/simulate.py` l. 128) | **yes** |
| A2 | decay between jumps | Euler, `y ← (1 − dt·β)·y` (`E/Simulate.f90` l. 49, 61; `AdditiveDrift = 0`, `E/SetParameters.f90` l. 21) | exact, `e^{−β dt}` (`simulate.py` l. 101) | yes (small) |
| A3 | jump arrival per step | Bernoulli with probability dt·λ (`E/Simulate.f90` l. 43, 55; `DiscreteDist1`, `E/Procedures.f90` l. 728–752) | probability 1 − e^{−λ dt} (`simulate.py` l. 102) | yes (small) |
| A4 | step length | dt = 0.25 quarter = **16 steps per year** (`E/Parameters.f90` l. 6) | 6 per quarter = **24 per year** (`simulate.py` l. 68) | yes |
| A5 | shock distribution | normal (`UseNormalDist = 1`, `E/SetParameters.f90` l. 18) | normal | no |

**B. Simulation design**

| # | element | KMV (file:line) | our v3 | differs? |
|---|---|---|---|---|
| B1 | initial condition | stationary normal, variance λσ²/(2β + λ(1 − ρ²)) (`E/Simulate.f90` l. 12–24, 31–34); at ρ = 0 the reset formula | **lifecycle**: everyone starts at z = 0 (`simulate.py` l. 104–105, `lifecycle=True`) | **yes** |
| B2 | burn-in | `Tburn = FLOOR(100/dt) + 1` = 401 steps = 100 quarters = **25 years** (`E/Parameters.f90` l. 7) | none | **yes** |
| B3 | panel kept | 20 quarters after the burn-in, `Tann = 5` annual observations (`E/Parameters.f90` l. 8–9; years built at `E/Simulate.f90` l. 74–78) | 36 years simulated, **31 used** (ages 25–55: `N_AGE_COLS = 31`, `scripts/pol_v3_estimate.py` l. 49–50) | **yes** |
| B4 | simulated workers | `nsim = 5000` (the comment says 50000) (`E/Parameters.f90` l. 5) | 50 000 in estimation, 100 000 out of sample (`pol_v3_estimate.py` l. 50–51) | yes |
| B5 | annual aggregation | sum of exp(z) over the 16 steps of a year, then log (`E/Simulate.f90` l. 71–79) | sum of exp(z) over the 24 steps, then log (`simulate.py` l. 119–128) | no (same principle) |
| B6 | random numbers | all drawn once before estimation (common random numbers), seed (7755, 7744), Intel Fortran `RANDOM_NUMBER` (`E/SetParameters.f90` l. 74–92) | common random numbers through a fixed numpy seed (1234) | no in principle; the streams cannot be reproduced |

**C. The eight moments**

| # | element | KMV (file:line) | our v3 | differs? |
|---|---|---|---|---|
| C1 | variance of log earnings | **one cross-section**: year 1 only, 1/n (`E/ComputeMoments.f90` l. 10–11) | mean over the 31 age columns of the cross-sectional variance (`simulate.py` l. 148) | **yes** |
| C2 | 1-year change | **only** year 2 − year 1, one change per worker (`E/ComputeMoments.f90` l. 40–43) | all consecutive pairs at ages 25–55, pooled (`simulate.py` l. 149) | **yes** |
| C3 | "5-year" change | year 5 − year 1 = **a 4-year gap** (`E/ComputeMoments.f90` l. 55–58) | a true 5-year gap, all pairs at ages 25–55, pooled (`simulate.py` l. 150) | **yes** |
| C4 | kurtosis | standardized 4th central moment with 1/n, μ₄/μ₂² (`E/ComputeMoments.f90` l. 41–48, 56–63) | the same (`simulate.py` l. 152–154) | no |
| C5 | fractions of small changes | share of \|Δ1\| strictly below 0.1 / 0.2 / 0.5, on the single year-1→2 change (`E/ComputeMoments.f90` l. 70–73) | the same thresholds and strict `<`, on all pooled Δ1 (`simulate.py` l. 162–164) | definition no; sample yes (C2) |
| C6 | targets | hard-coded US numbers 0.7, 0.23, 17.8, 0.46, **11.55**, 0.54, 0.71, 0.86 (`E/SetParameters.f90` l. 49–59; Table 2 prints 11.6) | `targets/pol_grid_v3.json` | n/a |

**D. Objective**

| # | element | KMV (file:line) | our v3 | differs? |
|---|---|---|---|---|
| D1 | residuals | relative, model/target − 1 (`E/dfovec.f90` l. 130, 135, 146, 151, 163, 173, 178, 183) | relative, (m − t)/t (`kmv_earnings/estimate.py` l. 101) | no |
| D2 | weights | **0.5 for kurt Δ1 and kurt Δ5**, 1 otherwise, applied as √w inside the residual (`E/dfovec.f90` l. 126, 145, 162) | 1 for all 8 (`estimate.py` l. 93–94, `weights=None`) | **yes** |
| D3 | norm | sum of squares F = Σ fᵢ² of the weighted residuals (`E/newuob-h.f` l. 344–349 and `f_value`, l. 622–629). The `objective fun` printed in `iterations.txt` is √(weighted mean of squares) (`E/dfovec.f90` l. 187–190), only printed, not minimized | Σ wₘ·relₘ² (`estimate.py` l. 102) | no (squares in both) |

**E. Optimizer**

| # | element | KMV (file:line) | our v3 | differs? |
|---|---|---|---|---|
| E1 | algorithm | DFLS: NEWUOA-H, a derivative-free trust-region least-squares method (Zhang, Conn & Scheinberg) (`E/newuoa-h.f` l. 1–14; called at `E/Estimation.f90` l. 110) | differential evolution + bounded Nelder-Mead polish (`estimate.py` l. 158, 167), and TikTak (`kmv_earnings/tiktak.py`) | **yes** |
| E2 | starting point | one start, `MakeGuess`: λ₁ 0.08, β₁ 0.7, σ₁ 1.6, λ₂ 0.007, β₂ 0.01, σ₂ 1.6 (`E/MakeGuess.f90` l. 9–19) | 5 DE seeds with KMV Table 3 injected into the population, + TikTak with 4096 Sobol points and the same warm start (`pol_v3_estimate.py` l. 63–64, 100, 106) | yes |
| E3 | restarts / convergence | `ndfls = 3` sequential restarts, each from the previous solution; `rhobeg = 5` and `rhoend = 1e-4` in the transformed space; `maxfun = 500·(n + 1)` = 3500 per restart (`E/SetParameters.f90` l. 61–63; `E/Estimation.f90` l. 101–112). The paper run used 441 evaluations in total (`E/earnings_estimation_paper/iterations.txt`) | DE `tol 1e-6`, `maxiter 200`, `popsize 15`; NM `xatol 1e-4`, `fatol 1e-8`, 800 iterations | yes |
| E4 | parameter space | logistic transforms: λ = 2·logistic(x), σ = 2·logistic(x), β = 1·logistic(x), i.e. **λ ∈ (0, 2), σ ∈ (0, 2), β ∈ (0, 1)** per quarter, for both components (`E/SetParameters.f90` l. 65–69; `E/dfovec.f90` l. 25, 46, 82, 89, 96, 103) | log box (`estimate.py` l. 33–39): λ₁ [0.01, 0.6], **β₁ [0.05, 8]**, σ [0.2, 3.5], λ₂ [0.001, 0.1], β₂ [1e-4, 0.2] | **yes**: β₁ ≤ 1 and σ ≤ 2 bind for Poland if v3's β₁ = 2.8 and σ₁ = 2.0 survive the other changes |
| E5 | choice of the estimate | the end point of the single run | the lowest out-of-sample f across starts (100k workers, another seed) | yes |
| E6 | standard errors | none (§29.2, §29.7) | bootstrap + asymptotic | yes (irrelevant here) |

**F. The Markov chain**

| # | element | KMV (file:line) | our v3 | differs? |
|---|---|---|---|---|
| F1 | size | **3 × 11** = 33 (`D/Parameters.f90` l. 5–6) | 5 × 15 = 75 (`kmv_earnings/discretize.py` l. 151) | yes |
| F2 | grid points | symmetric, end points ±width/2, power-spaced on each half with x^{1/k}, k = `gridpar` (`D/TransitionMatrix.f90` l. 36–53; `PowerSpacedGrid`, `D/Procedures.f90` l. 448–483) | bin edges power-spaced on ±3.5σ, points = conditional means of N(0, σ²) within each bin (`discretize.py` l. 78–91) | yes |
| F3 | how the grid is chosen | width and `gridpar` of each component (4 numbers) are **fitted by DFLS** (4 restarts, `rhobeg 5`) to the 8 moments of the simulated chain, against **the data targets** (`D/SetParameters.f90` l. 47–62; bounds and guesses `D/MakeGuess.f90` l. 8–20: width ∈ [σ, 10σ], start 2σ; `gridpar` ∈ [0.3, 1], start 0.7; `D/dfovec.f90` l. 59–126, same 0.5 kurtosis weights at l. 82, 99). Appendix D.1 says the grid matches "the combined continuous-state process"; the code matches the data targets | a fixed rule; no fitting | **yes** |
| F4 | jump landings | normal CDF over mid-point bins, centred at ρ·(origin) = 0, rows renormalized (`D/TransitionMatrix.f90` l. 75–123; ρ₁ = ρ₂ = 0 at `D/SetParameters.f90` l. 26, 31) | CDF over our bin edges, centred at the origin (additive) (`discretize.py` l. 112) | **yes** (reset vs additive) |
| F5 | drift | 2-point linear interpolation toward (1 − β)·y (`DriftPointsVersion = 2`, `D/SetParameters.f90` l. 13; `D/TransitionMatrix.f90` l. 125–173). The rate to the inner neighbour is β·y / (y − y_inner): **the same upwind scheme** as ours | upwind (`discretize.py` l. 104) | no |
| F6 | variance correction | none; the fitted width plays that role through the moments | analytic rescale to λσ²/2β per component (`discretize.py` l. 217), then `match_var_log` (l. 220–232) | yes |
| F7 | how the chain is simulated for the fit | Euler step `I + dt·Q`, dt = 0.25 (`D/Simulate.f90` l. 21–22); start from the chain's ergodic distribution (l. 56–57); 20 quarters = 5 years (`D/Parameters.f90` l. 9; `D/Simulate.f90` l. 68–71); nsim 50 000 (`D/Parameters.f90` l. 7); uniform CRN, seed (7755, 7744) (`D/Estimation.f90` l. 12–19); the same moments as C1–C5 (`D/ComputeMoments.f90`) | `expm(Q)` quarterly, lifecycle from the grid point nearest 0 | yes |
| F8 | exported object | the generator `ymarkov_combined` = Kronecker sum (`D/CombinedProcess.f90` l. 28–32), states sorted ascending (l. 39–55); also `ytrans_qu` = (I + dt·Q)⁴ (`D/Simulate.f90` l. 37–43) | `Q` (Kronecker sum) and `P = expm(Q)` | no for Q |
| F9 | parameters used | rounded to 4–5 digits (`D/SetParameters.f90` l. 23–31: 0.07987, 0.7616, 1.7352, 0.00656, 0.00939, 1.5259) | the estimate itself | minor |

**Things noticed in passing** (no effect on the results):
- `D/CombinedProcess.f90` l. 40 calls `MINLOC(ygrid_combined)` (the
  not-yet-filled global) instead of `lygrid_combined`. It returns index 1,
  which is the minimum anyway, so the sort is right.
- `D/TransitionMatrix.f90` l. 79–116 read the global `y1grid` / `y2grid`
  instead of the arguments `ly1grid` / `ly2grid`. They are the same arrays
  at the only call site (`D/dfovec.f90` l. 59).
- Appendix Table D.1 prints the discretized kurt Δ5 as 13.2;
  `D/earnings_discretization_paper/moments.txt` has 12.23. The other
  discretized entries agree with the file.
- The package ships two estimation runs. `E/earnings_estimation_paper/` is
  the paper (441 evaluations). `E/earnings_estimation_output/` is a later
  re-run (309 evaluations): λ₁ 0.07989, β₁ 0.7544, σ₁ 1.7341, λ₂ 0.006528,
  β₂ 0.009206, σ₂ 1.5273. The ReadMe (`R/ReadMe.txt`, "EARNINGS PROCESS
  ESTIMATION", step 1) attributes the gap to simulation error. It is KMV's
  own measure of run-to-run noise: β₁ −0.9%, β₂ −2.0%, λ₂ −0.5%.

**Summary.** The method differs in six substantive places: A1 (reset), B1–B3
(ergodic cross-section, 25-year burn-in, 5 years), C1–C3 (one cross-section,
4-year "Δ5"), D2 (kurtosis weight 0.5), E1–E4 (DFLS from one guess, β ≤ 1,
σ ≤ 2) and F1–F7 (a 3 × 11 grid fitted to the moments). Steps A2–A4 and B4 are
numerical details. Point 4 (§31.4) measures each one separately.

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
- **[Poland: DONE as v3, §26–§30, 2026-11-09; FRA/GER/USA deliberately not
  recomputed.]** **Robustness test for the level target (§24.3/§24.4): `std_resid` + model
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
- ~~**Chain size for step 2:** 75 states is the current default; may need to drop
  to 33 if the model solve is too slow.~~ **Settled (§28.3):** 5×15 solves in
  ≈ 1.1 s; keep it.
- **Jump specification: additive (ours, 2016 WP) vs reset (KMV AER, KMV
  code) — §29.5, updated 2026-11-09 (§29.7–29.8).** The package confirms
  reset in KMV's code. Our baseline stays additive (working-paper
  specification). The reset robustness run leaves the dispersion conclusion
  unchanged (1.43 vs 1.51) and halves z₂ persistence (16 vs 31 yr). Left for
  the supervisor: confirm the additive baseline.
- **GER/USA `β₂` pins use the superseded KMV anchor (§29.9).** They were
  pinned to var₂ = 0.9103, our additive formula on KMV's parameters. KMV's
  published chain has 0.647 for the persistent component. FRA/GER/USA are
  not to be recomputed (brief, 2026-11-08), so this is for the supervisor:
  re-pin to 0.647, or keep and state the convention.
- **v3_end puts β₂ at its lower bound (§26.8).** Robustness only, no chain.
  Report as "the end-of-window data favour even more persistence"; no action
  unless the supervisor wants v3_end as a HANK input.
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
