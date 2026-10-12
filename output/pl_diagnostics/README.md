# Poland diagnostics — Stage A (raw GRID data, before trusting the estimation)

Companion to `../../DECISIONS.md` Task 1 (2026-10-11). Built by
`scripts/pl_diagnostics.py --stage A`. All numbers below are from the raw
year-by-year GRID cells (`year_panel_*.csv` in this folder), **not** the
window-averaged targets in `targets/pol_grid.json` — the point of this stage
is to look at the data before it gets collapsed into eight averaged numbers.

Stage B (model fit, chain diagnostics, cross-country parameter plot) will be
added once the Poland estimation run (`output/pol/`) finishes — it was still
running when this file was written.

## What's in each file

- **`01_moments_by_year.png`** — 6 panels, PL vs USA/FRA/GER by year: var(log
  annual earnings), var(1y change), var(5y change), kurt(1y change),
  kurt(5y change), N obs (1y-change cell). Grey band = 2008–10.
- **`02_fractions_by_year.png`** — the three "fraction of small 1y changes"
  targets, same layout.
- **`03_nobs_level_by_year.png`** — sample size of the level-earnings cell,
  log scale, all four countries on one axis (checks for sample
  redefinitions / discontinuities).
- **`04_qq_shape.png`** — quantile-vs-normal plots (empirical quantile of the
  1y change vs the normal quantile at the same probability, with a
  same-variance Gaussian reference line) for Poland in 2002/2009/2016, and one
  year each for France and the USA. A straight line = Gaussian; an S-shaped
  curve bowing away from the line at the ends = leptokurtosis.
- **`year_panel_{usa,fra,ger,pol}.csv`** — the underlying per-year numbers.

## What looks fine

- **Sample size is rock-stable.** PL's `N obs` for both the level cell (~4.2–4.5M)
  and the 1y-change cell (~4.0–4.5M) moves smoothly year to year across the
  entire 2000–2020 window — no jumps, no discontinuities, no sign of a sample
  redefinition. Same for USA/FRA/GER. This is the check the task asked for
  specifically, and it passes cleanly.
- **Leptokurtosis is visible by eye in every country and every year checked**
  (`04_qq_shape.png`), including Poland in a pre-crisis year (2002), the
  crisis trough (2009), and a recent year (2016) — the same flat-middle,
  fat-tail S-shape KMV's jump-drift process is built to capture, not a
  crisis-specific artifact.
- **Poland's muted response to the 2008–10 crisis is economically sensible,
  not a red flag.** `var_d1` for PL does **not** spike in 2008 — it was
  already declining from a 2006 peak (0.306) and continues declining smoothly
  through 2008–2010 (0.283 → 0.255). Contrast USA, which jumps sharply in 2008
  (`var_d1` 0.34 → 0.38, `kurt_d1` drops 12.6 → 11.5 — a broad shock
  compressing the relative tail), and FRA/GER, which show smaller versions of
  the same pattern. Poland is the one GRID-covered EU economy that avoided a
  2009 recession ("zielona wyspa"), and the earnings-risk data agree with that
  history. Worth a sentence in the thesis as a sanity check, not a worry.

## What stands out and should get a sentence (not hidden)

- **A level-dispersion decline that runs through, not around, the crisis.**
  `var(log annual earnings)` for Poland rises 2002→2007 (0.496 → 0.566), then
  declines essentially monotonically through 2020 (→ 0.463) — a ~18% peak-to-
  trough move over 13 years. This is a slow secular trend, not a crisis break
  (it turns over a year *before* 2008). It means the single window-averaged
  `var_log_earns` target (0.508, full 2000–2020 window) sits roughly in the
  middle of a range the data moved through, rather than describing any single
  year well — worth flagging as a limitation of averaging over a 21-year
  window with a visible trend. GER moves in the *opposite* direction over the
  same decades (rising from ~0.55 to ~0.68) — an opposite secular pattern
  worth a footnote if GER and POL are being compared side by side.

- **A kurtosis and small-change-fraction regime shift that starts exactly at
  the crisis and lasts for years — the most notable feature in the Polish
  series.** `kurt_d1` for Poland troughs at 13.2–13.6 in 2006–2008, then jumps
  to 15.4 in 2009 and **stays elevated at 15–17 for seven years** (2009–2016)
  before declining back to 14.5 by 2019. `frac_d1_lt_10%` shows the mirror
  image: a dip around 2006–2008 (~0.40–0.43), then a rise to a 2013 peak of
  0.52, then decline. USA, FRA and GER all show a kurtosis *dip* in 2008
  (broad-based shock) that *reverts within 2–3 years* — Poland shows a dip
  that reverts quickly (by 2009) and then **overshoots upward for most of a
  decade**. That is a different, longer-lived pattern than any of the
  reference countries, coinciding with — but outlasting — the global
  financial crisis by far more than a typical business-cycle effect would.
  Two readings, neither confirmed here: (a) a genuine post-2009 change in how
  earnings risk is distributed in Poland (e.g. more "most people's earnings
  barely move, but a growing few see large jumps" — plausibly linked to
  structural change in the Polish labour market, EU-accession-era emigration
  of workers at the margin, or industry recomposition after 2009); or (b) a
  methodological break in how GRID constructs the residualized change series
  for Poland partway through the window (worth checking against the GRID
  codebook / documentation if a definitive answer matters for the write-up).
  **This is exactly the kind of thing the task asked to surface rather than
  average over** — the full-window target blends two visibly different
  sub-periods (2000–2008 vs 2009–2019). A pre/post-2009 split estimation would
  settle which reading is right; flagged as an open item, not done here.

- **Both patterns together mean the full-2000–2020 window average is a blend
  of two regimes**, not a stable single process over the period. This doesn't
  invalidate Task 1 (the instruction was explicitly to use the full window),
  but it is worth one sentence in the thesis, and a split-sample robustness
  check is cheap if the supervisor wants it.

## Nothing else jumped out

No single-year outliers, no implausible sign flips, no missing-data patterns
beyond the structurally expected ones (level variance absent for 2000–2001;
5-year-change statistics absent for 2016–2020, since a 5-year-ahead change
can't be computed within the sample; the 2020 row has no 1y/5y change
statistics at all, presumably the last available year at download time). None
of these affect the 2000–2020 average, since missing cells are skipped rather
than treated as zero.

## Stage B (pending)

Added once `output/pol/` has a finished estimate:
- 8-moment fit (data vs continuous model vs discretized chain)
- simulated 1y-change distribution overlaid on the empirical GRID percentiles
- chain stationary distribution (log and level), sample 36-year worker paths,
  var(log y) by age for all four countries
- 6-parameter cross-country comparison with KMV marked (bootstrap CIs added
  later if Task 3 finishes in time)
