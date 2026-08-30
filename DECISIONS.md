# Modelling decisions and their justifications

Running log of choices made in the project, with the reasoning behind each, so
they can be written up in the thesis rather than reconstructed from memory.

Format per entry: **what was decided → why → how it was verified → what it
implies for the write-up.**

---

## 0. Scope: the thesis target is Poland; USA/FRA/GER are method validation

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

**For the write-up.** State this openly and early. A reader who knows KMV will
notice the difference; better to have explained it than to be asked.

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
- **Residual kurtosis gap in validation** (14.5 simulated vs KMV's 16.5 at their
  own parameters), attributed to GKOS sample details not modelled here —
  principally the minimum-earnings threshold.
- **Fractions of small changes are interpolated**, not observed. GRID reports
  quantiles, not `P(|Δ| < c)`; each percentile column gives a point of the CDF
  and the fractions follow from `F(c) − F(−c)` by linear interpolation. This is
  interpolation (the thresholds lie inside the percentile range), never
  extrapolation, and is the intended use of GRID's quantile statistics.
  Accuracy depends on the fineness of the percentile grid, which is high here.

**For the write-up.** A short "limitations" subsection. Stating these is
stronger than leaving them to be discovered.

---

## 16. Step-2 interface notes (from reading the Auclert `annual-review` code)

> Step 2 has its own running log: **`../hank_step2/STEP2_LOG.md`**. It records
> the baseline reproduction (S1), the resolved frequency question (S2 —
> quarterly, confirmed six ways), and the one-asset-vs-two-asset assessment (S3
> — the one-asset model reproduces KMV's 20/80 direct–indirect split almost
> exactly). The notes below are the pre-step-2 reading of the code and still
> stand.

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

## Open decisions (not yet settled)

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
- **USA `β₂` is on the high side too.** The unrestricted USA estimate implies an
  ergodic sd of the persistent component of ≈ 1.22 (vs KMV's 0.95), so the
  exported USA chain has ergodic var(log e) ≈ 1.8 — nearly double KMV's. Not
  implausible enough to force the §8b treatment, but if the step-2 USA steady
  state looks off, applying the same ergodic-var anchor to USA (for consistency
  with GER) is the fallback.
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
- **One-asset vs two-asset model in step 2:** ~~open~~ **largely resolved**
  (STEP2_LOG §S3). The one-asset ARS model (β-heterogeneity) reproduces KMV's
  headline 20/80 direct–indirect split for the monetary shock almost exactly
  (19.7% direct / 80.3% indirect, first-year, vs KMV's 19/80). The *composition*
  of the indirect effects differs (asset-return channel +22% vs KMV's −2%,
  transfer channel smaller) but that is the non-household blocks, not one- vs
  two-asset. Proceed one-asset. Two-asset only needed if the thesis wants the
  cross-sectional MPC distribution (KMV Fig 5–6), the wealthy-HtM mechanism as
  such, or the portfolio channel. Supervisor's call.
- **Poland:** earnings targets (Polish panel data / EU-SILC) and wealth moments
  (NBP HFS / ECB HFCS) still pending.
