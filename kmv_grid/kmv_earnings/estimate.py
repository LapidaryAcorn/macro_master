"""
SMM estimation of the KMV earnings process.

Objective: sum of squared *relative* deviations of the 8 simulated moments
from the data targets (all targets are strictly positive, so relative
deviations put the variance-type and kurtosis-type moments on a common scale).

Optimiser: scipy differential_evolution (global) over log-parameters,
optionally polished with Nelder-Mead. Common random numbers inside the
simulator keep the objective stable across evaluations.
"""

from __future__ import annotations

import json
import numpy as np
from scipy.optimize import differential_evolution, minimize

from .simulate import MOMENT_ORDER, PARAM_ORDER, model_moments

# bounds on quarterly rates (log-scale search inside these)
DEFAULT_BOUNDS = {
    "lambda1": (0.01, 0.6),
    "beta1":   (0.05, 4.0),
    "sigma1":  (0.20, 3.5),
    "lambda2": (0.001, 0.10),
    "beta2":   (0.001, 0.20),
    "sigma2":  (0.20, 3.5),
}


def _vec_to_params(x: np.ndarray) -> dict:
    return {name: float(np.exp(v)) for name, v in zip(PARAM_ORDER, x)}


def make_objective(targets: dict, weights: dict | None = None, sim_kwargs: dict | None = None):
    sim_kwargs = sim_kwargs or {}
    w = np.array([1.0 if weights is None else weights.get(m, 1.0) for m in MOMENT_ORDER])
    t = np.array([targets[m] for m in MOMENT_ORDER])

    def objective(x: np.ndarray) -> float:
        params = _vec_to_params(x)
        m = model_moments(params, **sim_kwargs)
        mv = np.array([m[k] for k in MOMENT_ORDER])
        rel = (mv - t) / t
        return float(np.sum(w * rel**2))

    return objective


def estimate(
    targets: dict,
    bounds: dict | None = None,
    weights: dict | None = None,
    sim_kwargs: dict | None = None,
    maxiter: int = 60,
    popsize: int = 12,
    seed: int = 7,
    polish_nm: bool = True,
    workers: int = 1,
    disp: bool = True,
    x0: dict | None = None,
):
    """
    Run global SMM estimation. Returns (params_hat, result_object).

    For a serious run use e.g. sim_kwargs=dict(n_workers=50_000, n_years_keep=10)
    and maxiter >= 100. For quick pipeline tests, shrink n_workers/maxiter.
    `x0`: optional dict of starting parameters injected into the initial population
    (e.g. the KMV US estimates as a warm start for FR/DE).
    """
    bounds = bounds or DEFAULT_BOUNDS
    log_bounds = [tuple(np.log(bounds[p])) for p in PARAM_ORDER]
    obj = make_objective(targets, weights=weights, sim_kwargs=sim_kwargs)

    init = "latinhypercube"
    if x0 is not None:
        rng = np.random.default_rng(seed)
        pop = np.array([
            [rng.uniform(lo, hi) for lo, hi in log_bounds]
            for _ in range(popsize * len(PARAM_ORDER))
        ])
        pop[0] = [np.log(x0[p]) for p in PARAM_ORDER]
        init = pop

    res = differential_evolution(
        obj, log_bounds, maxiter=maxiter, popsize=popsize, seed=seed,
        init=init, tol=1e-6, mutation=(0.4, 1.0), recombination=0.7,
        polish=False, disp=disp, workers=workers, updating="deferred" if workers != 1 else "immediate",
    )

    x_best = res.x
    if polish_nm:
        nm = minimize(obj, x_best, method="Nelder-Mead",
                      options={"xatol": 1e-4, "fatol": 1e-8, "maxiter": 400})
        if nm.fun < res.fun:
            x_best = nm.x

    params_hat = _vec_to_params(x_best)
    return params_hat, res


def save_params(params: dict, path: str) -> None:
    with open(path, "w") as f:
        json.dump(params, f, indent=2)


def load_params(path: str) -> dict:
    with open(path) as f:
        return json.load(f)
