# KMV 2018 earnings process: estimation on GRID data (Table 2 replication)

Pipeline: GRID moments -> SMM estimation of the KMV (2018, AER) two-component
jump-drift process -> "Earnings Process Estimation Fit" table (Data / Model
Estimated / Model Discretized) -> export of the 33-state quarterly Markov
chain for the HANK / sequence-space model (step 2).

## Layout
- `kmv_earnings/simulate.py`   - process simulation + the 8 moments
- `kmv_earnings/estimate.py`   - SMM (differential evolution + Nelder-Mead polish)
- `kmv_earnings/discretize.py` - 3x11=33-state Markov chain, moments, txt export
- `kmv_earnings/grid_loader.py`- building the 8 targets from GRID csv files
- `kmv_earnings/run.py`        - CLI
- `targets/`                   - target moment JSONs (US from KMV Table 2 included)

## Quick start
```bash
pip install numpy scipy pandas

# 1. sanity check of the simulator at KMV's published Table 3 parameters
python -m kmv_earnings.run validate

# 2. estimation (US benchmark first; then FR/DE with GRID targets)
python -m kmv_earnings.run estimate --targets targets/us_kmv.json \
    --out output/us --warm-start --workers 4

# quick pipeline test only (small sim, few iterations):
#   add --quick

# 3. rebuild table + chain export from saved parameters
python -m kmv_earnings.run table --targets targets/us_kmv.json \
    --params output/us/params.json --out output/us
```
Outputs per run: `table2.csv`, `table2.tex` (the Table-2-style table),
`params.json`, and `income_process_{grid,zgrid,P,pi}.txt` (quarterly 33-state
chain, mean earnings normalised to 1) - the input for step 2.

## Model and design decisions (validated against the published paper)
Process (KMV eq. 30-31): log z = z1 + z2, dz_j = -beta_j z_j dt + eps_j dN_j,
eps_j ~ N(0, sigma_j^2), Poisson arrivals lambda_j. Quarterly rates.
Three details matter a lot and were verified by matching KMV's fitted moments
at their published Table 3 estimates (lambda=[0.080, 0.007],
beta=[0.761, 0.009], sigma=[1.74, 1.53]):
1. Jumps are ADDITIVE (not resets of the component).
2. Lifecycle panel: workers enter at z=0 and are followed 36 years
   (~ages 25-60, the GKOS SSA sample). With lambda2=0.007 the persistent
   shock arrives ~once per career, so the ergodic distribution would
   overstate observed inequality (var(log y) ~0.9 vs 0.70).
3. Simulate at high frequency (default 6 steps/quarter) and time-aggregate
   the earnings FLOW exp(z) to annual; moments are on log annual earnings.
With these, simulation at KMV's parameters gives
0.69 / 0.23 / 0.50 / 14.6 / 10.8 / 0.55 / 0.66 / 0.84
vs KMV's model column
0.70 / 0.23 / 0.46 / 16.5 / 12.1 / 0.56 / 0.67 / 0.85.
The residual kurtosis gap plausibly reflects GKOS sample details (minimum
earnings threshold, entry/exit) not modelled here.

The "Model Discretized" column uses a self-contained scheme (conditional-mean
grid points, upwind drift, additive-jump generator, expm). It matches KMV's
discretized variances/fractions well but understates kurtosis; for an exact
match, swap in the discretization routine from the KMV replication package
(appendix D.1) - the rest of the pipeline is unaffected.

## Getting targets from GRID (https://www.grid-database.org)
Download for a country (e.g. DE, FR): (i) log earnings LEVELS statistics,
(ii) 1-year and (iii) 5-year residualized log-change statistics, males,
prime age, averaged over years. Then either:
- fill a JSON by hand (copy `targets/us_kmv.json`), or
- use `grid_loader.targets_from_grid_csvs(levels_csv, d1_csv, d5_csv)`.

Notes:
- GRID does NOT report "frac |change| < 10/20/50%" directly; the loader
  interpolates them from the percentile grid of the 1y-change distribution.
- Use the standard (moment-based) kurtosis, not Crow-Siddiqui, for
  comparability with KMV's 17.8 / 11.6.
- CAVEAT: the loader's column names (`COLUMN_MAP`) were written without live
  access to the GRID portal - check them against the actual csv headers and
  the GRID data dictionary after downloading.
- Comparability with KMV's US targets: males, ~25-55(60), residualized
  changes, minimum-earnings threshold. Match the GRID cells accordingly.
- Test #1 (as planned): build US targets from GRID and check you recover
  KMV's Table 2/3; then FR and DE; then PL when NBP/supervisor data arrive.

## Estimation practicalities
- 6 parameters, 8 moments; objective = sum of squared relative deviations.
- Common random numbers (fixed seed, fixed RNG-call structure) keep the
  objective deterministic across parameter values.
- Serious run: defaults (50k workers, 36y, 6 steps/q), maxiter >= 100,
  --workers N for parallel DE. `--warm-start` seeds the population with the
  KMV US estimates (sensible for FR/DE).
- Identification intuition: var/kurt at 1y vs 5y horizons separate the
  transitory and persistent components; fractions of small changes pin down
  arrival rates; var(log y) pins overall scale.

## TikTak vs differential evolution (US targets, equal budgets ~950 evals)
See `output/comparison.txt` for the full run. Summary: in-sample TikTak
f=0.0276 vs DE f=0.0345; re-evaluated out-of-sample on a large simulation
(150k workers, 6 steps/q, fresh seed) they are a statistical tie
(0.048 vs 0.047), both slightly below KMV's own Table 3 parameters under
this simulator (0.056). CAVEAT: at these short budgets the objective is
flat in some directions - DE's solution has a near-random-walk persistent
component (lambda2=0.047, beta2~0, sigma2=0.40) that fits the 8 moments as
well as a KMV-like configuration but has very different economic content.
TikTak's solution stays close to KMV. For thesis results: run longer
(n_sobol >= 512, more locals / maxiter >= 100), warm-start from KMV, and
check that the persistent-component parameters are stable across methods.
CLI: `python -m kmv_earnings.run estimate --method tiktak ...`

## Real GRID export format (verified July 2026)
`data/Stats_*.csv`: one wide file, rows = country x year (e.g. USA/FRA/GER),
columns: std_log_inc, std/kurt of residualized 1y and 5y log changes, and a
fine percentile grid p1...p99_99 of the 1y change. Use
`grid_loader.targets_from_grid_stats_csv(path, country)` -> 8-moment dict.
IMPORTANT: check the `gender` column of your export. "All genders" is NOT
comparable with KMV's male-only GKOS targets (it inflates var(log y) and
var(dlog y) and lowers the small-change fractions). Re-download selecting
males only for the US validation test; e.g. an All-genders US export gave
var_log 0.93 vs KMV's 0.70.
