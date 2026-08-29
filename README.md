# KMV (2018) earnings process — estimation on GRID data

Step 1 of a master's thesis applying the Kaplan–Moll–Violante (2018, AER) HANK
framework cross-country. This repository estimates the KMV two-component
jump-drift earnings process on harmonized GRID data and exports a finite-state
Markov chain for the HANK / sequence-space model (step 2).

Pipeline: **GRID statistics → 8 target moments → SMM estimation → Table-2-style
fit table + exported quarterly Markov chain.**

---

## Layout

| Path | What it does |
|---|---|
| `kmv_earnings/simulate.py` | The income process itself + the 8 moments |
| `kmv_earnings/estimate.py` | SMM objective + differential evolution (+ Nelder-Mead polish) |
| `kmv_earnings/tiktak.py` | TikTak global optimiser (Arnoud–Guvenen–Kleineberg) |
| `kmv_earnings/discretize.py` | Markov-chain discretization, diagnostics, txt export |
| `kmv_earnings/grid_loader.py` | Builds the 8 targets from GRID `Stats_*.csv` |
| `kmv_earnings/run.py` | CLI (`validate` / `estimate` / `table`) |
| `KMV_step1_test.ipynb` | Step-by-step notebook covering the whole pipeline |
| `data/` | GRID export |
| `targets/` | Target-moment JSONs per country |

## Quick start

```bash
pip install -r requirements.txt          # numpy, scipy, pandas (+ jinja2 for .tex)

# 1. sanity check: simulate at KMV's published Table 3 parameters
python -m kmv_earnings.run validate

# 2. estimate on a country's GRID targets
python -m kmv_earnings.run estimate --targets targets/fra_grid.json \
    --out output/fra --warm-start
#    --quick          small sim + few iterations (pipeline test only)
#    --method tiktak  use TikTak instead of differential evolution
#    --workers N      parallel DE

# 3. rebuild table + chain export from saved parameters
python -m kmv_earnings.run table --targets targets/fra_grid.json \
    --params output/fra/params.json --out output/fra
```

Outputs per run: `table2.csv`, `table2.tex`, `params.json`, and
`income_process_{grid,zgrid,P,pi}.txt` — the quarterly Markov chain (mean
earnings normalised to 1) that step 2 consumes.

---

## The model

KMV eq. 30–31: log earnings `z = z1 + z2`, each component

```
dz_j = -beta_j * z_j dt + eps_j dN_j,    eps_j ~ N(0, sigma_j^2),  N_j ~ Poisson(lambda_j)
```

Component 1 is transitory (frequent, fast-decaying), component 2 persistent
(rare, near-permanent). Six parameters, quarterly rates.

**Three design details matter a lot** and were each verified by reproducing
KMV's fitted moments at their published Table 3 estimates
(`lambda=[0.080, 0.007]`, `beta=[0.761, 0.009]`, `sigma=[1.74, 1.53]`):

1. **Jumps are additive**, not resets of the component.
2. **Lifecycle panel**: workers enter at `z = 0` and are followed 36 years
   (≈ ages 25–60, the GKOS SSA sample). With `lambda2 = 0.007` the persistent
   shock arrives roughly once per career, so the ergodic distribution would
   overstate observed inequality (var(log y) ≈ 0.9 vs 0.70).
3. **High-frequency simulation with time aggregation**: simulate at 6 steps per
   quarter and aggregate the earnings *flow* `exp(z)` to annual; all moments are
   computed on log annual earnings.

With these, simulating at KMV's parameters gives

```
          var_log  var_d1  var_d5  kurt_d1  kurt_d5   f<10%  f<20%  f<50%
ours         0.69    0.23    0.50    14.5     10.6     0.55   0.66   0.84
KMV model    0.70    0.23    0.46    16.5     12.1     0.56   0.67   0.85
```

The residual kurtosis gap plausibly reflects GKOS sample details (minimum
earnings threshold, entry/exit) not modelled here. `python -m kmv_earnings.run
validate` reproduces this table.

---

## Discretization

The HANK model needs a finite state space, so each component is put on a grid
(bin edges power-spaced, grid points = conditional means of `N(0, sigma^2)`
within each bin), the drift is an upwind finite-difference scheme, additive
jumps integrate the normal density over bins, and the quarterly transition
matrix is `expm(Q)`. The combined chain is the Kronecker product of the two
independent components.

**Two corrections were needed** (both diagnosed against the continuous
benchmark, see the `discretize.py` docstring for the full write-up):

*Variance bias.* The upwind scheme is numerically diffusive and inflates the
chain's ergodic variance by ~19% (1.27 vs the exact `lambda*sigma^2/(2*beta)`
total of 1.07). The bias is systematic — it does **not** shrink as the grid is
refined (33 → 465 states leaves it unchanged). Fixed analytically by rescaling
each component's grid by `sqrt(v_theory / v_chain)`; kurtosis is scale-invariant,
so this corrects the variances without touching the tails. `diagnose()` reports
chain-vs-theory variances per component.

*Grid too coarse for the tails.* Kurtosis needs grid points far enough out to
represent rare large jumps. KMV's 3×11 = 33 states cannot carry it; the default
is now **5×15 = 75 states**.

Measured against the continuous simulation (kurt_d1 = 14.6):

| grid | states | var_log | var_d1 | kurt_d1 | kurt_d5 |
|---|---|---|---|---|---|
| (3, 11) | 33 | 0.63 | 0.19 | 10.9 | 9.5 |
| **(5, 15)** | **75** | 0.64 | 0.20 | **13.9** | 10.2 |
| (9, 21) | 189 | 0.64 | 0.20 | 14.7 | 10.6 |
| (11, 25) | 275 | 0.64 | 0.20 | 14.8 | 10.7 |

Past ~(9, 21) nothing is gained. Pass `n1`/`n2` explicitly if step 2 needs a
smaller chain (33 states still benefits from the variance fix: kurt 11.9 vs 8.1
before) or can afford a larger one.

**Known residual gaps, honestly stated.** (a) Variances land ~7% below the
continuous process, because the analytic correction targets the *ergodic*
variance while moments come from a 36-year lifecycle panel — pass
`match_var_log=True` to close this by simulation. (b) `frac_d1_lt_10` comes out
too *high* (~0.65 vs 0.55): the chain often does not change state at all in a
quarter, piling mass at near-zero annual changes. This is intrinsic to
discretization and no grid refinement removes it. (c) For an exact match to
KMV's published discretized column, substitute the routine from their
replication package (appendix D.1); nothing else depends on how `P` is built.

---

## Targets from GRID

GRID export format (verified against a real download): one wide
`Stats_*.csv`, rows = country × year, with `std_log_inc`, `std`/`kurt` of
residualized 1-year and 5-year log changes, and a fine percentile grid
`p1 … p99_99` of the 1-year change.

```python
from kmv_earnings.grid_loader import targets_from_grid_stats_csv, save_targets
t = targets_from_grid_stats_csv("data/Stats_20260720145608.csv", "FRA")
save_targets(t, "targets/fra_grid.json")
```

**What to download:** males, prime age (25–55), residualized changes. You need
log-earnings *levels* statistics, plus 1-year and 5-year change statistics —
and the **percentiles of the 1-year change**, which are not optional.

**Fractions of small changes are interpolated, not read off.** GRID does not
report `P(|Δ| < 10/20/50%)`; it reports quantiles. Since the quantile function
is the inverse CDF, each percentile column gives a known point of the CDF, and
`P(|Δ| < c) = F(c) − F(−c)` is obtained by linear interpolation between
neighbouring percentiles (`fractions_from_percentiles`). All three thresholds
are interpolated — none coincides with a reported percentile. This is
interpolation, never extrapolation, and it is the intended use of GRID's
quantile statistics. Accuracy depends on how fine the percentile grid is, so
download it as fine as offered.

Use the standard moment-based kurtosis, **not** Crow–Siddiqui, for
comparability with KMV.

### Current data: coverage and the KMV comparison

Male, 25–55. Year coverage differs by country: **FRA 1991–2016, GER 2001–2016,
USA 1998–2019.**

|  | KMV (GKOS) | USA | FRA | GER |
|---|---|---|---|---|
| var(log annual earns) | 0.70 | 0.945 | 0.488 | 0.645 |
| var 1yr change | 0.23 | 0.323 | 0.206 | 0.148 |
| var 5yr change | 0.46 | 0.618 | 0.336 | 0.285 |
| kurt 1yr change | 17.8 | 12.87 | 15.86 | 17.80 |
| kurt 5yr change | 11.6 | 8.83 | 11.68 | 11.21 |
| frac 1yr < 10% | 0.54 | 0.387 | 0.528 | 0.611 |
| frac 1yr < 20% | 0.71 | 0.587 | 0.725 | 0.788 |
| frac 1yr < 50% | 0.86 | 0.811 | 0.864 | 0.894 |

**GRID-USA does not match KMV's published US column, and this is expected.**
KMV target the GKOS SSA extract over ~1978–2013; GRID-USA here is 1998–2019 with
GRID's own harmonized sample construction. US earnings dispersion rose over
these decades, so a later window gives higher variances and correspondingly
lower fractions of small changes — and that is exactly the pattern: every
dispersion moment is inflated and every small-change fraction is depressed, by
amounts that shrink as the threshold widens (−0.15 at 10%, −0.05 at 50%). All
eight moments tell one consistent "more volatile sample" story, which is the
signature of a sample-definition difference rather than a bug. The likely second
contributor is GKOS's minimum-earnings threshold, which trims volatile
low-earners.

The validity standard for the thesis is therefore **consistency of method across
countries**, not replication of KMV's specific US extract — and this should be
stated explicitly in the text.

**On a common year window.** Restricting all countries to their overlap
(2001–2016) barely moves anything (USA var_log 0.945 → 0.956, FRA 0.488 → 0.476,
GER unchanged). The trends within these windows are mild enough that coverage
differences are not driving the cross-country comparison. Pass
`year_range=(2001, 2016)` to `targets_from_grid_stats_csv` if you want the
strictly comparable version; either way, state the choice in the thesis.

---

## Estimation

Six parameters, eight moments. The objective is the sum of squared *relative*
deviations, `Σ ((simulated − target)/target)^2`, so that variances (~0.2) and
kurtoses (~15) carry comparable weight. The search runs in log-parameters, and
the simulator uses **common random numbers** (fixed seed, fixed RNG-call
structure), which makes the objective a deterministic, smooth-ish function of
the parameters instead of a noisy one.

Bounds keep the two components in their intended roles — component 1 frequent
and fast-decaying, component 2 rare and near-permanent — which prevents the
optimiser from relabelling or collapsing them. Check that no estimate sits *at*
a bound; if one does, widen it and re-run.

Two global optimisers are available:

- **Differential evolution** (`--method de`, default): evolves a population via
  difference-vector mutation; step size self-adapts to the population spread.
  Followed by a Nelder-Mead polish.
- **TikTak** (`--method tiktak`): Sobol space-filling scan, then a sequence of
  local searches started from convex combinations of each seed and the incumbent
  best, with the weight on the incumbent rising along the sequence.

**Run both and compare — this is an identification check, not a contest.** The
persistent component (`lambda2`, `beta2`, `sigma2`) is weakly identified: with
`lambda2 ≈ 0.007` a 36-year panel contains few persistent shocks, so several
configurations fit the 8 moments almost equally well. If the two optimisers land
in the same place, the estimate is trustworthy; if they diverge, first run longer
budgets, and if they *still* diverge, the problem is identification rather than
search — the fix is an additional identifying moment (e.g. a longer-horizon
autocovariance, which GRID reports), not more computation.

Observed on USA (equal budgets, out-of-sample scoring on a large fresh-seed
panel): the two agree closely, including on the persistent component
(`lambda2` 0.010 vs 0.012, `beta2` 0.008 vs 0.007), objectives 0.085 vs 0.093.
On earlier short-budget runs DE occasionally found a near-random-walk persistent
component that fit equally well but meant something economically different —
which is precisely why the check exists.

For thesis numbers use the full defaults (not `--quick`), `maxiter >= 100` or
`--n-sobol 512`, and `--warm-start`.

---

## Status and next steps

**Done.** Validated simulator; verified GRID loader; targets for USA/FRA/GER;
SMM estimation with two global optimisers; corrected discretization; Table-2
output and Markov-chain export; test notebook.

**Next.**
1. Full-budget estimations for FRA and GER, with the DE-vs-TikTak agreement
   check for each country.
2. Regenerate the exported chains — the files currently in `output/` were built
   with the old 33-state default and predate the discretization fix.
3. Decide and document the year window.
4. **Step 2**: `shade-econ/annual-review` (Auclert et al. sequence-space
   Jacobian). Start from the one-asset model and replace the Rouwenhorst
   discretization of a lognormal AR(1) with loading `income_process_*.txt`.
   Watch: mean-1 normalisation (already handled), frequency (their model is
   quarterly, as is the chain — an annual variant would need `P^4`), and that
   75 states instead of the usual 7–11 will slow the solve. Sanity path:
   reproduce their steady state with the original process, swap in the chain,
   compare the wealth distribution and MPCs, then run the monetary IRFs.
5. Poland: earnings targets from Polish panel data / EU-SILC, wealth moments
   from NBP HFS / ECB HFCS.

## References

Kaplan, Moll & Violante (2018, AER), *Monetary Policy According to HANK* — the
core paper. At the earnings process they cite: Guvenen, Karahan, Ozkan & Song
(2015/2021) for the leptokurtic moments and the symmetry assumption; Kaplan &
Violante (2014) for the two-asset structure; Schmidt (2015) for compound-Poisson
income dynamics; Karahan & Ozkan (2013) on shock persistence. Also relevant:
Guvenen, Ozkan & Song (2014, JPE) on countercyclical income risk; Arnoud,
Guvenen & Kleineberg (2019) for TikTak; Guvenen (2009, RED) on identification of
persistent components; McFadden (1989) / Pakes & Pollard (1989) for SMM.
