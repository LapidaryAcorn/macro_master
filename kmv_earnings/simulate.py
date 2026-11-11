"""
Simulation of the KMV (2018, AER) earnings process and computation of the
8 target moments from Table 2 ("Earnings Process Estimation Fit").

Model (KMV eq. 30-31), quarterly units:
    log earnings  z_it = z1_it + z2_it
    dz_j = -beta_j * z_j dt + eps_j dN_j,   eps_j ~ N(0, sigma_j^2)
i.e. each component decays toward zero at rate beta_j and receives ADDITIVE
jumps eps_j arriving at Poisson rate lambda_j.
Component 1 = transitory (high lambda, high beta), 2 = persistent.

Sample design (validated against KMV's fitted moments):
  - lifecycle panel: workers enter at z = 0 and are followed for 36 years
    (~ages 25-60, the GKOS SSA sample). With lambda2 = 0.007 the persistent
    shock arrives about once per career, so the ergodic distribution would
    overstate observed inequality.
  - simulate at high frequency (default 6 steps per quarter) and time-
    aggregate the earnings FLOW exp(z) to annual earnings; all moments are
    computed on log annual earnings.

At KMV's Table 3 estimates this reproduces their model column closely:
sim: 0.69 0.23 0.50 14.6 10.8 0.55 0.66 0.84
KMV: 0.70 0.23 0.46 16.5 12.1 0.56 0.67 0.85
(residual kurtosis gap likely reflects sample details of GKOS, e.g. the
minimum-earnings threshold, not implemented here).
"""

from __future__ import annotations

import numpy as np

MOMENT_ORDER = [
    "var_log_earns", "var_d1", "var_d5", "kurt_d1", "kurt_d5",
    "frac_d1_lt_10", "frac_d1_lt_20", "frac_d1_lt_50",
]

MOMENT_LABELS = {
    "var_log_earns": "Variance: annual log earns",
    "var_d1": "Variance: 1yr change",
    "var_d5": "Variance: 5yr change",
    "kurt_d1": "Kurtosis: 1yr change",
    "kurt_d5": "Kurtosis: 5yr change",
    "frac_d1_lt_10": "Frac 1yr change < 10%",
    "frac_d1_lt_20": "Frac 1yr change < 20%",
    "frac_d1_lt_50": "Frac 1yr change < 50%",
}

PARAM_ORDER = ["lambda1", "beta1", "sigma1", "lambda2", "beta2", "sigma2"]

# KMV (2018) Table 3, quarterly rates (verified against the published paper).
KMV_TABLE3_PARAMS = {
    "lambda1": 0.080, "beta1": 0.761, "sigma1": 1.74,
    "lambda2": 0.007, "beta2": 0.009, "sigma2": 1.53,
}

# KMV Table 2, "Data" column (US, SSA via Guvenen-Karahan-Ozkan-Song).
KMV_US_DATA_TARGETS = {
    "var_log_earns": 0.70, "var_d1": 0.23, "var_d5": 0.46,
    "kurt_d1": 17.8, "kurt_d5": 11.6,
    "frac_d1_lt_10": 0.54, "frac_d1_lt_20": 0.71, "frac_d1_lt_50": 0.86,
}


def simulate_log_annual_earnings(
    params: dict,
    n_workers: int = 50_000,
    n_years_keep: int = 36,
    steps_per_quarter: int = 6,
    lifecycle: bool = True,
    burn_in_years: int = 0,
    seed: int = 1234,
    deterministic_arrival: bool = False,
    reset_jumps: bool = False,
    kmv_euler: bool = False,
    exact_reset_burn: bool = False,
) -> np.ndarray:
    """
    Return log annual earnings, shape (n_workers, n_years_keep).

    `kmv_euler`: the per-step rules of KMV's EarningsEstimation/Simulate.f90
    (l. 43, 49, 55, 61): a jump arrives with probability lambda*dt (not
    1 - exp(-lambda*dt)) and between jumps z decays by the Euler factor
    (1 - beta*dt) (not exp(-beta*dt)), dt = 1/steps_per_quarter quarters. With
    steps_per_quarter=4 this is KMV's 16 steps per year (DECISIONS Sec 31).
    Default False: unchanged, bit for bit (tests/test_kmv_method.py).

    Fixed seed + fixed RNG-call structure => common random numbers across
    parameter values (the SMM objective is smooth-ish and deterministic).

    `deterministic_arrival`: KMV's own footnote nests a discrete-time AR(1) in
    logs as the limit of their jump-drift process where shocks arrive every
    period with probability 1 (lambda -> infinity) rather than at a Poisson
    rate -- i.e. z_j(t+dt) = z_j(t) e^{-beta dt} + eps every sub-step, no
    Bernoulli draw. Setting this True reproduces that nested special case
    (Task 2, DECISIONS.md 2026-10-11); lambda_j is then irrelevant for
    component j and should be fixed to an arbitrary placeholder.

    `reset_jumps`: on arrival the component is REPLACED by eps, z_j -> eps,
    instead of z_j -> z_j e^{-beta dt} + eps. This is the process in KMV's
    replication code (rho = 0 in their Simulate.f90 / TransitionMatrix.f90)
    and the AER text, eq. 33-34 (DECISIONS Sec 29.5, 29.7). Stationary
    variance lam*sig^2/(lam + 2*beta). Default False = additive, bit-for-bit
    unchanged (same RNG calls; tests/test_reset_jumps.py).

    `exact_reset_burn` (needs reset_jumps, burn_in_years > 0): draw the state
    at the end of the burn-in from its exact distribution instead of stepping
    through it. Under resets with per-step jump probability p and decay d the
    state after B steps is z0*d^B if no jump occurred (probability (1-p)^B),
    else eps*d^K with K = steps since the last jump, P(K=k) prop. to
    p(1-p)^k, k < B. One uniform and one normal per worker and component,
    so common random numbers are kept; the distribution of every simulated
    path is the same as with step-by-step burn-in (tests/test_kmv_method.py),
    the random streams are not. About 6x faster in KMV's design (400 of 480
    steps are burn-in). Default False: unchanged.
    """
    lam = np.array([params["lambda1"], params["lambda2"]]) / steps_per_quarter
    beta = np.array([params["beta1"], params["beta2"]]) / steps_per_quarter
    sig = np.array([params["sigma1"], params["sigma2"]])

    rng = np.random.default_rng(seed)
    decay = np.exp(-beta)
    p_jump = np.ones(2) if deterministic_arrival else 1.0 - np.exp(-lam)
    if kmv_euler:
        decay = 1.0 - beta
        p_jump = np.ones(2) if deterministic_arrival else lam

    z = np.zeros((n_workers, 2))
    if not lifecycle:
        # stationary variance of the jump-drift process: lam*sig^2/(2*beta)
        # (additive) or lam*sig^2/(lam + 2*beta) (reset)
        for j in range(2):
            if reset_jumps:
                sd = sig[j] * np.sqrt(lam[j] / (lam[j] + 2 * beta[j]))
            else:
                sd = sig[j] * np.sqrt((lam[j] / (2 * beta[j])) if beta[j] > 0 else 0.0)
            z[:, j] = rng.standard_normal(n_workers) * sd

    spy = 4 * steps_per_quarter                     # steps per year
    total_s = spy * (burn_in_years + n_years_keep)
    burn_s = spy * burn_in_years

    if exact_reset_burn and burn_s > 0:
        if not reset_jumps or deterministic_arrival:
            raise ValueError("exact_reset_burn needs reset_jumps and Poisson arrival")
        B = burn_s
        for j in range(2):
            u = rng.random(n_workers)
            eps = rng.standard_normal(n_workers) * sig[j]
            p = float(p_jump[j])
            if p <= 0.0:
                z[:, j] = z[:, j] * decay[j] ** B
                continue
            l1p = np.log1p(-p)
            q = np.exp(B * l1p)                       # P(no jump in B steps)
            jumped = u >= q
            v = (u - q) / (1.0 - q)                   # uniform on [0, 1) given a jump
            k = np.ceil(np.log1p(-v * (1.0 - q)) / l1p) - 1.0
            k = np.clip(np.where(jumped, k, 0.0), 0.0, B - 1)
            z[:, j] = np.where(jumped, eps * decay[j] ** k, z[:, j] * decay[j] ** B)
        total_s -= burn_s
        burn_s = 0

    logy = np.empty((n_workers, n_years_keep))
    annual_flow = np.zeros(n_workers)
    year_idx = 0

    for s in range(total_s):
        for j in range(2):
            u = rng.random(n_workers)
            eps = rng.standard_normal(n_workers) * sig[j]
            zd = z[:, j] * decay[j]
            z[:, j] = np.where(u < p_jump[j], eps if reset_jumps else zd + eps, zd)
        if s >= burn_s:
            annual_flow += np.exp(z[:, 0] + z[:, 1])
            if (s - burn_s) % spy == spy - 1:
                logy[:, year_idx] = np.log(annual_flow)
                annual_flow[:] = 0.0
                year_idx += 1
    return logy


def moments_from_panel(logy: np.ndarray, n_age_cols: int | None = None,
                       d5_gap: int = 5, cross_section: bool = False) -> dict:
    """8 target moments from a (workers x years-of-age) panel.

    `n_age_cols`: if set, use only the first `n_age_cols` columns (ages), e.g.
    31 = ages 25-55 of the 36-year (25-60) panel, matching GRID's age window
    (DECISIONS Sec 26). 1y/5y changes are then taken within that window only.
    None (default) uses every column, exactly as before.

    `d5_gap`: the gap, in years, of the "5-year" change (default 5). KMV's
    code uses year 5 - year 1, i.e. 4 (EarningsEstimation/ComputeMoments.f90
    l. 55-58; DECISIONS Sec 31).
    `cross_section`: KMV's moments from ONE cross-section (ComputeMoments.f90
    l. 10-11, 40-43, 55-58, 70-73): the variance of column 0 only, one 1-year
    change (col 1 - col 0) and one "5-year" change (col d5_gap - col 0) per
    worker. Default False: pooled over all columns, unchanged.
    """
    if n_age_cols is not None:
        logy = logy[:, :n_age_cols]
    g = int(d5_gap)
    if cross_section:
        var_log = float(np.var(logy[:, 0]))
        d1 = logy[:, 1] - logy[:, 0]
        d5 = logy[:, g] - logy[:, 0]
    else:
        var_log = float(np.mean(np.var(logy, axis=0)))
        d1 = (logy[:, 1:] - logy[:, :-1]).ravel()
        d5 = (logy[:, g:] - logy[:, :-g]).ravel()

    def kurt(x):
        x = x - x.mean()
        return float(np.mean(x**4) / np.mean(x**2) ** 2)

    return {
        "var_log_earns": var_log,
        "var_d1": float(np.var(d1)),
        "var_d5": float(np.var(d5)),
        "kurt_d1": kurt(d1),
        "kurt_d5": kurt(d5),
        "frac_d1_lt_10": float(np.mean(np.abs(d1) < 0.10)),
        "frac_d1_lt_20": float(np.mean(np.abs(d1) < 0.20)),
        "frac_d1_lt_50": float(np.mean(np.abs(d1) < 0.50)),
    }


def model_moments(params: dict, n_age_cols: int | None = None, d5_gap: int = 5,
                  cross_section: bool = False, **sim_kwargs) -> dict:
    """`n_age_cols`, `d5_gap`, `cross_section` as in moments_from_panel. With
    a lifecycle panel the first k simulated years do not depend on how many
    years follow (same RNG calls, in the same order), so only
    min(n_years_keep, n_age_cols) years are simulated -- identical moments,
    less work (tests/test_age_window.py)."""
    if n_age_cols is not None and sim_kwargs.get("lifecycle", True):
        sim_kwargs = dict(sim_kwargs)
        sim_kwargs["n_years_keep"] = min(sim_kwargs.get("n_years_keep", 36), n_age_cols)
    return moments_from_panel(simulate_log_annual_earnings(params, **sim_kwargs),
                              n_age_cols=n_age_cols, d5_gap=d5_gap,
                              cross_section=cross_section)
