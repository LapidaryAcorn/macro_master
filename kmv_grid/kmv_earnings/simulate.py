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
) -> np.ndarray:
    """
    Return log annual earnings, shape (n_workers, n_years_keep).

    Fixed seed + fixed RNG-call structure => common random numbers across
    parameter values (the SMM objective is smooth-ish and deterministic).
    """
    lam = np.array([params["lambda1"], params["lambda2"]]) / steps_per_quarter
    beta = np.array([params["beta1"], params["beta2"]]) / steps_per_quarter
    sig = np.array([params["sigma1"], params["sigma2"]])

    rng = np.random.default_rng(seed)
    decay = np.exp(-beta)
    p_jump = 1.0 - np.exp(-lam)

    z = np.zeros((n_workers, 2))
    if not lifecycle:
        # stationary variance of the jump-drift process: lam*sig^2/(2*beta)
        for j in range(2):
            sd = sig[j] * np.sqrt((lam[j] / (2 * beta[j])) if beta[j] > 0 else 0.0)
            z[:, j] = rng.standard_normal(n_workers) * sd

    spy = 4 * steps_per_quarter                     # steps per year
    total_s = spy * (burn_in_years + n_years_keep)
    burn_s = spy * burn_in_years

    logy = np.empty((n_workers, n_years_keep))
    annual_flow = np.zeros(n_workers)
    year_idx = 0

    for s in range(total_s):
        for j in range(2):
            u = rng.random(n_workers)
            eps = rng.standard_normal(n_workers) * sig[j]
            zd = z[:, j] * decay[j]
            z[:, j] = np.where(u < p_jump[j], zd + eps, zd)
        if s >= burn_s:
            annual_flow += np.exp(z[:, 0] + z[:, 1])
            if (s - burn_s) % spy == spy - 1:
                logy[:, year_idx] = np.log(annual_flow)
                annual_flow[:] = 0.0
                year_idx += 1
    return logy


def moments_from_panel(logy: np.ndarray) -> dict:
    var_log = float(np.mean(np.var(logy, axis=0)))
    d1 = (logy[:, 1:] - logy[:, :-1]).ravel()
    d5 = (logy[:, 5:] - logy[:, :-5]).ravel()

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


def model_moments(params: dict, **sim_kwargs) -> dict:
    return moments_from_panel(simulate_log_annual_earnings(params, **sim_kwargs))
