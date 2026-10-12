# The estimated earnings processes: what they are, how well they fit, and how much to trust them

This note describes the KMV (2018) jump-drift earnings process estimated for
each of the four countries in this project — the United States, France,
Germany, and Poland — and states plainly how much confidence each parameter
deserves. It is written to be read start to finish, not consulted as a log.
The decisions behind each choice, with full derivations and dead ends, are
in `DECISIONS.md`; this note draws on that record (principally §§8–23) but
does not reproduce its reasoning.

## 1. The model and how it is fit, in brief

Log annual earnings are the sum of two independent components, each a
continuous-time jump-drift process: `dz_j = -β_j z_j dt + ε_j dN_j`, with
`ε_j ~ N(0, σ_j²)` arriving at Poisson rate `λ_j`. Component 1 is
transitory — frequent, fast-decaying; component 2 is persistent — rare,
slow-decaying. Six parameters per country: `λ₁, β₁, σ₁, λ₂, β₂, σ₂`, all
quarterly rates.

Each country's eight target moments — variance of log annual earnings,
variance and kurtosis of 1-year and 5-year log changes, and the fractions of
1-year changes below 10%/20%/50% — come from GRID (the Global Repository of
Income Dynamics), for male workers aged 25–55. The parameters are estimated
by simulated method of moments: minimize the sum of squared *relative*
deviations between simulated and target moments, on a 36-year lifecycle
panel (workers start at zero and are followed through a career, not drawn
from the process's stationary distribution — the two differ sharply when
persistent shocks arrive roughly once per working life). Two global
optimizers are run for every country — differential evolution (DE) and
TikTak (Sobol global phase + sequential Nelder-Mead local phase) — and DE is
the one actually relied on: TikTak's local phase stalls on this objective's
stair-stepped surface (common random numbers make the simulated moments
change in discrete jumps as parameters cross grid-cell boundaries), a
behavior confirmed directly rather than assumed (§9, §21).

Three further exercises test how much the resulting numbers can be trusted,
beyond the fit itself:

- **A bootstrap** over GRID's own year-to-year statistics (resample years
  with replacement, recompute the moments, re-estimate) gives standard
  errors — but only for year-to-year *sampling* variability in the reported
  statistics, and only *conditional on the point estimate's own basin*
  (§19, §21).
- **A parametric recovery exercise** takes each country's own estimate as
  ground truth, simulates a fresh synthetic panel from it with an
  independent seed, and re-estimates blind, to see whether the procedure
  can find its way back (§21) — a sharper identification test than checking
  the two optimizers agree, because it does not require already knowing
  where to look.
- **Nested simpler specifications** (one component only; deterministic
  rather than Poisson arrival, i.e. KMV's own AR(1) limiting case; a pure
  random walk persistent component) show which of the eight moments break
  when pieces of the jump-drift structure are removed (§22), and — for
  Poland specifically — whether that cost in *fit* corresponds to anything
  that matters once the chain is fed into the step-2 HANK model (§23).

## 2. The parameters, economically

| quantity | USA | FRA | GER | POL |
|---|---:|---:|---:|---:|
| Transitory half-life (years) | 0.059 | 0.041 | 0.033 | 0.061 |
| Persistent half-life (years) | 12.73 | 33.63 | 26.42 | 10.25 |
| Mean time between transitory shocks (years) | 1.16 | 1.31 | 1.58 | 1.54 |
| Mean time between persistent shocks (years) | 23.58 | 35.47 | 24.43 | 28.49 |
| σ, transitory | 1.941 | 1.910 | 1.791 | 1.979 |
| σ, persistent | 1.529 | 1.092 | 1.080 | 1.401 |
| Ergodic sd, transitory component | 0.373 | 0.286 | 0.218 | 0.334 |
| Ergodic sd, persistent component | 0.954 | 0.903 | 0.954 | 0.714 |
| Ergodic sd, combined | 1.024 | 0.948 | 0.979 | 0.788 |

(`tables/10_economic_params.csv` / `.tex`; USA and GER's persistent-component
rows reflect the ergodic-variance anchor they are *calibrated* to, not a
free estimate — see §4.1 and §4.3 below.)

Two things stand out immediately, independent of any later caveat. First,
the transitory component decays in **weeks, not quarters**, for every
country — half-lives of 12–22 days. A worker's short-run earnings wiggles
are gone almost as fast as they arrive; this is also why a quarterly
discretized chain systematically understates the 1-year change variance by
about 30% (§8d) — the transitory shock mostly decays within a single
quarter, faster than a quarterly Markov chain can resolve. Second, the
persistent component arrives roughly **once per working life** everywhere
(23–35 years between shocks) and, once it arrives, is slow enough that it is
effectively permanent over a career even where it is not literally a random
walk. This is why the 36-year lifecycle panel used for estimation and the
process's *ergodic* distribution are not interchangeable (§2) and why the
persistent component's decay rate is hard to pin down from this kind of
panel in the first place (§8b, §21, §23) — a worker observed for one career
has lived through only a handful of these shocks, so the panel carries very
little information about how fast they decay.

## 3. Fit quality

All four countries fit the eight target moments well on the continuous
model — largest relative error under 5% almost everywhere, **with one
specific, known exception**: USA and Germany, the two countries whose
persistent-component ergodic variance is *pinned* to KMV's anchor rather
than freely estimated (§4.1, §4.3), undershoot `var_log_earns` by about
**15%** on the continuous model (USA −16.0%, GER −14.6%; France and Poland,
both freely estimated, are within 0.6% of data on this moment). This is a
direct, expected consequence of the pin, not a separate fit problem — it is
the fit cost of choosing a believable long-run ergodic dispersion over
whatever value would best match the 36-year lifecycle-panel moment, and it
should be reported alongside the pin itself rather than left implicit in
the full fit tables (full tables: `output/{usa,fra,ger,pol}/table2.csv`;
Stage-B fit plots: `output/pl_diagnostics/05_model_fit_*.png`).

The discretized 75-state chain used for step 2 carries further known,
documented gaps that are a property of discretizing a quarterly chain, not
of the fit itself: the 1-year change variance comes in roughly 25–45% low
(worst for GER and FRA, whose transitory components decay fastest), and the
"fraction of tiny 1-year changes" comes in too high, because in many
quarters the chain simply does not move at all (§8d, §15). A further,
separately documented discretization gap — the 1-year change *kurtosis*
coming in substantially *higher* in the discretized chain than in the
continuous model (+18% to +41% depending on country) — is detailed in §15.
None of these four gaps is specific to one country.

## 4. Country by country

### 4.1 United States

The free, unrestricted 6-parameter estimate (`output/usa_free_beta2/`)
implies an ergodic standard deviation for the persistent component of
**1.22** — well above KMV's own 0.95 and above every other country here.
Confronting an infinite-horizon HANK model with that much ergodic earnings
risk breaks its calibration outright (STEP2_LOG §S4), so the **canonical**
USA estimate (`output/usa/`) pins the persistent component's ergodic
variance to KMV's own anchor (0.9103) via `ErgodicVarConstraint` and
estimates the remaining five parameters freely (§8b, applied to USA in the
open-decisions history). The table above reports the canonical, pinned
version.

**Recovery (§21):** on the *free* estimate, `β₂` recovers to within 1–4% of
truth across 4 independent synthetic draws (2 seeds × 2 panel sizes), with
no degradation at the larger panel size — a clean, well-identified result.
**Nested specs (§22):** forcing `β₂=0` costs USA 3.6× in OOS fit
(`f`: 0.0058 → 0.0212) — a real but modest cost. **Bootstrap SE**
(`output/bootstrap/`, 200 reps): `λ₁=0.214±0.004`, `β₁=2.97±0.08`,
`σ₁=1.94±0.01`, `λ₂=0.0104±0.0003`, `σ₂=1.56±0.02`; `β₂=0.0140` is reported
without an SE because it is calibrated, not estimated, in the canonical run.
No line scan was run specifically for USA (§9's line scan covers FRA and
GER, the two countries where the two optimizers visibly diverged on
`β₁`; USA's two optimizers agreed to begin with).

The GRID-USA sample itself is also a known outlier in *shape*, separate
from the ergodic-variance issue: thinner-tailed (`kurt Δ1y` = 12.8 against
15–18 for the other three countries) and with a slower-decaying transitory
component (`β₁` the lowest of the four at the canonical estimate). STEP2_LOG
§S4–S5 document that this shape, not just the variance, is what makes the
GRID-USA chain fail to support the KMV transmission mechanism in a one-asset
HANK model even after the variance is pinned — the step-2 US validation
therefore runs on KMV's own published chain, with GRID-USA reported as a
sensitivity whose failure is itself informative (§12).

### 4.2 France

France is estimated with the full 6 free parameters throughout — no pin was
ever needed; its free ergodic sd (0.90) already sits close to KMV's own
anchor. **Line scan (§9):** the two optimizers landed at different `β₁`
(DE 4.27, TikTak-recheck 2.18), but a scan along the straight line between
them in log-parameter space is strictly monotonically decreasing from
TikTak's point to DE's — a single valley, DE simply further down it, not two
competing basins. Extending that same line **25–100% past the DE point**
(§9, 2026-10-12) shows the objective rising sharply — more than doubling by
25% further, over 40-fold by 100% further — confirming DE's point is the
floor of the valley, not a point still on a descending slope. **Recovery
(§21):** `β₂` recovers to within 3.5–4.4% across 4 draws, tight and stable —
well identified, with a small (and, per the Monte Carlo noise check in §9,
partly simulation-noise-sized) bias that does not shrink with sample size.
**Nested specs (§22):** `β₂=0` costs France 8.3× in OOS fit
(0.0013 → 0.0108) — the largest proportional cost of the three
non-GER countries, though still a small absolute number. **Bootstrap SE**:
`λ₁=0.190±0.002`, `β₁=4.37±0.14`, `σ₁=1.91±0.01`, `λ₂=0.0069±0.0002`,
`β₂=0.0048±0.0004`, `σ₂=1.094±0.003` — note `β₂`'s own SE is about 8% of its
point estimate, the largest relative SE of any free parameter reported here,
consistent with it being the hardest parameter to pin down even where it is
not literally unidentified.

### 4.3 Germany

Germany is the one country in this project where a parameter is
**calibrated rather than estimated**, and the evidence for that choice is
now about as direct as this kind of evidence gets. The free 6-parameter
estimate (`output/ger_free_beta2/`) puts `β₂ ≈ 0.0006` — a half-life near
250 years, for practical purposes a random walk — with an implied ergodic
sd of **2.69**, an order of magnitude more dispersed in the long run than
any other country and not a credible description of any labour market. The
canonical estimate instead pins the persistent component's ergodic variance
to KMV's own 0.9103 (§8b) and estimates the other five parameters freely;
that is the row reported in §2's table.

**Why pin rather than estimate.** Three independent pieces of evidence, not
one:
1. *An objective sweep* (§8b) already showed only a weak interior optimum
   near `β₂ ≈ 0.0007`, i.e. an almost flat direction.
2. *The line scan* (§9): DE and TikTak disagree sharply on `β₁` (5.34 vs
   2.27); the line between them is monotone (bar one stair-step within
   simulation noise) from TikTak to DE, and — extended 25–100% past DE
   (2026-10-12) — rises clearly and monotonically past it, quadrupling by
   100% further. So the canonical point is a genuine local minimum, not an
   open slope, *along this one direction*.
3. *Recovery* (§21) is the sharpest evidence: simulating noise-free
   synthetic data from the free estimate and re-estimating blind, with a
   full serious budget and a neutral (not truth-informed) start, **DE lands
   at a different point in all 4 of 4 independent draws** — `β₁ ≈ 1.84`
   (vs the true 5.30) and `β₂ ≈ 0.0001` (vs the true 0.0006), every time,
   regardless of seed or synthetic panel size. The consistency *is* the
   finding: this is a repeatable feature of the objective surface, not
   optimizer noise. Checking the objective value at the true point versus
   this alternate point on the same synthetic data shows the true point is
   genuinely better (`f=0.00024` vs `0.00085`) — so it is a flat valley the
   optimizer cannot climb out of with a serious-but-finite budget, not two
   literally tied optima. Strikingly, the alternate point DE converges to
   sits in the same low-`β₁` neighbourhood where TikTak's own Nelder-Mead
   search stalled on GER's **real** data (§9) — a different algorithm, a
   different (synthetic vs real) dataset, landing in the same place.
   Checked directly whether this casts doubt on the canonical free estimate
   itself: on GER's real data, the canonical estimate still beats this
   alternate point by 2.6× in objective terms — a clearer margin than on
   synthetic data, so the canonical estimate stands, but the margin is not
   enormous.

**What this is not evidence of:** `β₁`, `λ₁`, `σ₁`, `λ₂`, `σ₂` are all well
recovered in the same exercise (§21) — the identification problem is
specific to `β₂` (and, in the free/unrestricted specification only, drags
`β₁`/`λ₁` down with it into the same alternate basin; the canonical, pinned
specification does not have this problem because `β₂` is not free in it).

**Nested specs (§22):** forcing `β₂=0` costs GER essentially nothing in fit
— `f: 0.0025 → 0.0028`, a 1.1× ratio that, per the Monte Carlo noise
measured directly in §9 (~23% of `f`'s own value at this scale), is
indistinguishable from zero cost. This is exactly consistent with the free
estimate already sitting near `β₂=0`. **Bootstrap SE** (five free
parameters; `β₂=0.0066` reported without an SE, being calibrated):
`λ₁=0.156±0.002`, `β₁=5.52±0.13`, `σ₁=1.80±0.01`, `λ₂=0.01029±0.00009`,
`σ₂=1.083±0.002` — note these SEs describe variability *within the
canonical, pinned basin only* (§19); they say nothing about the risk,
demonstrated directly in §21, that an unconstrained search lands somewhere
else entirely.

### 4.4 Poland

Poland is the newest estimate and gets the most space here because it is
the country whose data arrived last and whose robustness needed the most
direct checking before anything is built on it.

**Data and estimation.** GRID's Polish export uses different column names
than the FRA/GER/USA file (`std_perm`/`_resid_dNyr` instead of
`std_log_inc`/`_res_Nyr_log_chg`); the loader now auto-detects both schemas,
and the mapping was checked semantically, not just for successful parsing —
a hand-computed reproduction of six moments directly from the raw CSV,
using nobody else's helper code, matched the loader's output exactly
(`DECISIONS.md` §18). Targets use Poland's full native 2000–2020 window
(not the 2001–2016 window shared by the other three countries — a
deliberate choice, not an oversight; see the window-robustness check below).
The full 6-parameter estimate is free throughout, like France. DE and
TikTak **agree** on the persistent component to within 1.3% — no red flag
from that check, for what it is worth (§9 already cautions this check alone
is not conclusive for any country).

**Where Poland sits, shape-wise.** Its transitory half-life (0.061 years)
and persistent half-life (10.25 years, the *shortest* of the four
countries) put it closer to the USA on some dimensions and closer to
France/Germany on others — it is not a clean match to either reference
country, which matters because the pre-existing plan in this project was to
use "does Poland look USA-shaped or FRA/GER-shaped" as a trigger for
whether a one-asset HANK model would suffice (see the resolved
open-decisions item below).

**A real post-2009 regime shift in the raw data, checked and found not to
change any downstream decision.** The raw year-by-year GRID statistics show
`kurt Δ1y` jumping from ~13.5 to 15–17 starting exactly in 2009 and staying
elevated for seven years, at roughly unchanged variance
(`output/pl_diagnostics/README.md`) — a genuine shape change, not noise,
and not matched by anything in USA/FRA/GER's own 2008–10 response. Splitting
the sample at 2009 and re-estimating both halves (`scripts/pol_regime_split.py`)
shows the *parameters* move only modestly despite the sharp moment shift —
mainly a somewhat less active persistent component after 2009 (lower
`λ₂`/`β₂`/`σ₂`) — and, critically, **the two-asset trigger check below gives
the identical answer on both sub-windows and the full window** (§20). So the
2000–2020 window choice, while it does blend two visibly different
sub-periods, turns out not to be consequential for anything this project
currently does with the estimate.

**Two-asset trigger: resolved, and the rule of thumb would have gotten it
wrong.** The pre-existing plan was that a low `β₁` (≲3) or thin tail would
signal Poland needs a two-asset HANK model rather than the one-asset +
death structure that works for France and Germany. Poland's `β₁ = 2.85`
sits just under that threshold. Running the actual check
(`hank_step2/swap_death.py pol`) rather than applying the rule mechanically
shows Poland calibrates **cleanly** — converged, `MPC = 0.2000` exactly,
`frac HtM = 0.226`, `dβ = 0.125` (healthy, nowhere near collapsing), direct/
indirect decomposition 22.6%/77.4% — matching France and Germany, nothing
like the USA's degenerate collapse (§20). The rule of thumb pointed the
wrong way; the actual calibration run did not.

**Recovery (§21): mostly excellent, with one real but diagnosable
exception.** 3 of 4 independent draws recover `β₂` to within 1%. The
4th (the largest synthetic panel, 300,000 workers) collapses to the same
qualitative alternate basin seen systematically in Germany — but checking
the objective at the true point versus the recovered point on that draw's
synthetic data shows an 83-fold gap (`f=0.00011` vs `0.00947`), far larger
than Germany's 3.5-fold gap on an analogous check. That asymmetry is the
tell: this is DE failing to converge on one unlucky run, not a robust
competing optimum — Poland's identification margin is real, just not as
comfortable as France's or the USA's, and repeating the estimate from more
than one starting point (already standard practice here via the DE/TikTak
check) remains worth keeping specifically for Poland.

**Nested specs (§22), and why the small fit cost is misleading here.**
Forcing `β₂=0` costs Poland only 3.5× in OOS fit (0.0037 → 0.0130) — on its
face the second-cheapest restriction after Germany's. But Poland's `V3`
estimate sits at or against two parameter bounds (`λ₂ = 0.09999636`
against its 0.10 ceiling, `σ₂ = 0.213` against its 0.20 floor) and moves
`β₁` from 2.85 down to 0.50 — a more than fivefold change in the
*transitory* component. The restricted estimate is not "the same process
with `β₂` switched off"; it is a structurally different process that
happens to reproduce the eight moments about as well. The 3.5× fit-cost
number understates how much this restriction actually costs Poland.

**The ergodic-variance profile (§23): the keystone result for any
PL-vs-USA comparison, and the reason this note exists.** `β₂` itself never
enters the step-2 HANK model; what does, via the discretized chain's
ergodic distribution, is `var2 = λ₂σ₂²/(2β₂)`. Pinning `var2` at four
values and re-estimating the rest at full budget
(`scripts/pol_beta2_profile.py`), then running the full `swap_death.py`
calibration on each resulting chain, shows the calibration **fails to
converge at both ends of the tested range** — not only the low end
(`var2=0.30`, where `β`-heterogeneity hits its bound and MPC overshoots to
0.32), but also at `var2=0.9103`, **the exact ergodic-variance anchor
already imposed on the USA and Germany** for the unrelated reason that
their own free estimates do not support a clean calibration (MPC overshoots
to 0.33 there too). Only `var2` in roughly `[0.4, 0.7]` — bracketing
Poland's own free estimate of 0.509 — calibrates cleanly, and that range is
itself internally indistinguishable from noise (§23a: `f` at `var2=0.40`
and `var2=0.70` differ by only 1.5%, well inside the ~23% Monte Carlo noise
band measured directly in §9). The asymmetric discipline is clean and
moment-specific on both sides: `var_log_earns` pulls `var2` up from below
(worst moment at both `var2=0.30` and `0.40`), `var_Δ5y` pulls it down from
above (worst moment at both `var2=0.70` and `0.9103`).

**The consequence, stated plainly because it belongs in the main text, not
an appendix:** Poland's one-asset-plus-death calibration is clean precisely
*because* it runs on Poland's own data-implied ergodic dispersion. Forcing
it onto the same anchor already used for the USA — the natural move if the
thesis wants a PL-vs-USA comparison built on a common calibration target —
breaks it. Any comparison across the two countries built that way is
confounded: part of whatever difference it finds will reflect Poland
running on its own dispersion while the USA runs on an imposed one, not a
genuine difference in process shape alone.

**Bootstrap SE** (200 reps, all 6 parameters free): `λ₁=0.159±0.003`,
`β₁=2.850±0.055` (95% CI [2.74, 2.96], tightly centered on the point
estimate 2.846 — direct confirmation that an earlier validation scare,
traced to a stale warm-start file from a race condition with a concurrently
running job, was a methodological artifact, not a real problem, §19),
`σ₁=1.982±0.011`, `λ₂=0.00882±0.00035`, `β₂=0.0169±0.0012`,
`σ₂=1.411±0.016`.

## 5. Cross-country synthesis

**On shape, Poland does not cleanly match either reference pattern.** Its
transitory decay (`β₁=2.85`) sits closer to the USA (2.75–2.97 depending on
restriction) than to France or Germany (4.3–5.5); its persistent-component
ergodic dispersion (sd 0.71) is the *lowest* of the four countries, not the
highest, and its kurtosis (data 15.2) sits between the USA's thin tail
(12.8) and France/Germany's fatter ones (15.2–17.8). This mixed profile is
exactly why the two-asset trigger question could not be settled by a rule
of thumb and needed the direct check in §20.

**On identification, the four countries split three ways, not two.** USA
and France recover `β₂` well (1–4% in synthetic recovery) and pay only a
modest fit cost to restrict it. Germany's `β₂` is genuinely unidentified —
recovery fails systematically, in every draw, landing on a known alternate
basin, and the restriction costs nothing in fit because the data do not
distinguish the restricted and unrestricted specifications at all. Poland
sits closer to USA/France (good recovery, small-sample identification) but
is demonstrably not as robust as either — one draw in four falls into the
same qualitative trap Germany falls into every time, even though the
objective-gap evidence says that one draw is an optimizer failure rather
than a second real basin.

**On what matters downstream, the fit-based identification story and the
HANK-based one diverge for Poland specifically.** Sections §21/§22 would,
taken alone, suggest `β₂`'s weak fit-discipline is a second-order nuisance
for every country except Germany. §23 shows this is not quite right for
Poland: the SMM objective is close to indifferent across a `var2` range of
roughly `[0.4, 0.7]`, but the HANK calibration is not indifferent to where
in — or outside — that range the chain actually sits, and in particular is
not indifferent to being pushed onto the same anchor used for the USA. That
distinction — "cheap in the moment-fitting sense" versus "consequential
once fed into the model that is actually used" — is the single most
important methodological point this round of work adds.

iMPC sequences (`hank_step2/results/impc_comparison.png`,
`impc_profile.png`) and the capital-gains-response sequences
(`mcap`) for all six valid chains plus the four Poland ergodic-variance
profile points are available as the primary cross-country comparison
object recommended by the literature check in `DECISIONS.md` §17
(Auclert–Rognlie–Straub's sufficient-statistic result, read carefully: it
applies exactly to fixed-real-rate fiscal shocks, and only via a
generalization — requiring log utility and no cross-type portfolio
heterogeneity, both satisfied by this project's one-asset model — to the
interest-rate shocks this project's step 2 actually studies).

## 6. What this note does not cover

Wealth moments for Poland (NBP HFS / ECB HFCS) are still pending — the
`beta` calibration target for a Polish-specific HANK run is a separate,
unresolved thesis-design question (`DECISIONS.md`, open-decisions list).
The channel-level (as opposed to aggregate 20/80) decomposition match to
KMV remains imperfect for every country, for reasons traced to the
non-household GE blocks rather than the income process (STEP2_LOG §S3).
Nested specifications were only run for the four canonical chains; the same
exercise has not been repeated on either of Poland's sub-window estimates
or on the ergodic-variance profile points.
