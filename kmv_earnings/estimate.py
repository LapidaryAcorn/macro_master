"""
SMM estimation of the KMV earnings process.

Objective: sum of squared *relative* deviations of the 8 simulated moments
from the data targets (all targets are strictly positive, so relative
deviations put the variance-type and kurtosis-type moments on a common scale).

Optimiser: scipy differential_evolution (global) over log-parameters,
polished with a *bounded* Nelder-Mead. Common random numbers inside the
simulator keep the objective stable across evaluations.

Restricting the parameter space
-------------------------------
`estimate` / `tiktak` accept `fixed` (hold a parameter at a constant) and
`derive` (compute a parameter from the free ones each evaluation). Both shrink
the search vector. `derive` values must be *picklable* callables (so parallel
`workers` still work on Windows/spawn) - use `ErgodicVarConstraint` or another
module-level class, not a lambda. Example (GER persistent component pinned so
its ergodic variance equals KMV's):

    derive = {"beta2": ErgodicVarConstraint("2", target_var=0.9103)}
"""

from __future__ import annotations

import json
import numpy as np
from scipy.optimize import differential_evolution, minimize

from .simulate import MOMENT_ORDER, PARAM_ORDER, model_moments

# bounds on quarterly rates (log-scale search inside these)
DEFAULT_BOUNDS = {
    "lambda1": (0.01, 0.6),
    "beta1":   (0.05, 8.0),    # widened: FRA/GER transitory decay wants ~4-4.5
    "sigma1":  (0.20, 3.5),
    "lambda2": (0.001, 0.10),
    "beta2":   (0.0001, 0.20),  # lowered: near-permanent components sit well below 1e-3
    "sigma2":  (0.20, 3.5),
}


class ErgodicVarConstraint:
    """Picklable `derive` callable: return beta_j such that component j's
    ergodic variance  lambda_j * sigma_j^2 / (2 beta_j)  equals `target_var`.
    Keeps a weakly identified near-permanent component's long-run dispersion at
    a chosen plausibility level instead of letting it run to implausible values.
    """

    def __init__(self, comp: str = "2", target_var: float = 0.9103):
        self.lam = f"lambda{comp}"
        self.sig = f"sigma{comp}"
        self.target_var = float(target_var)

    def __call__(self, p: dict) -> float:
        return p[self.lam] * p[self.sig] ** 2 / (2.0 * self.target_var)


def _free_params(fixed: dict | None, derive: dict | None) -> list[str]:
    locked = set(fixed or {}) | set(derive or {})
    return [p for p in PARAM_ORDER if p not in locked]


def _build_params(x: np.ndarray, free: list[str],
                  fixed: dict | None = None, derive: dict | None = None) -> dict:
    params = {name: float(np.exp(v)) for name, v in zip(free, x)}
    if fixed:
        params.update({k: float(v) for k, v in fixed.items()})
    if derive:
        for k, fn in derive.items():
            params[k] = float(fn(params))
    return params


def _vec_to_params(x: np.ndarray) -> dict:
    """Back-compat: full 6-vector in PARAM_ORDER -> params dict."""
    return {name: float(np.exp(v)) for name, v in zip(PARAM_ORDER, x)}


class _Objective:
    """SMM objective as a top-level (picklable) callable, so scipy's
    differential_evolution can farm it out to worker processes (needed on
    Windows / spawn, where a closure cannot be pickled)."""

    def __init__(self, targets: dict, weights: dict | None = None,
                 sim_kwargs: dict | None = None, free: list[str] | None = None,
                 fixed: dict | None = None, derive: dict | None = None):
        self.sim_kwargs = sim_kwargs or {}
        self.fixed = fixed or {}
        self.derive = derive or {}
        self.free = (list(free) if free is not None
                     else _free_params(self.fixed, self.derive))
        self.w = np.array([1.0 if weights is None else weights.get(m, 1.0)
                           for m in MOMENT_ORDER])
        self.t = np.array([targets[m] for m in MOMENT_ORDER])

    def __call__(self, x: np.ndarray) -> float:
        params = _build_params(x, self.free, self.fixed, self.derive)
        m = model_moments(params, **self.sim_kwargs)
        mv = np.array([m[k] for k in MOMENT_ORDER])
        rel = (mv - self.t) / self.t
        return float(np.sum(self.w * rel**2))


def make_objective(targets: dict, weights: dict | None = None,
                   sim_kwargs: dict | None = None, free: list[str] | None = None,
                   fixed: dict | None = None, derive: dict | None = None):
    return _Objective(targets, weights=weights, sim_kwargs=sim_kwargs,
                      free=free, fixed=fixed, derive=derive)


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
    fixed: dict | None = None,
    derive: dict | None = None,
):
    """
    Run global SMM estimation. Returns (params_hat, result_object).

    For a serious run use e.g. sim_kwargs=dict(n_workers=50_000, n_years_keep=36)
    and maxiter >= 100. For quick pipeline tests, shrink n_workers/maxiter.
    `x0`: optional dict of starting parameters injected into the initial
    population (e.g. the KMV US estimates as a warm start for FR/DE).
    `fixed` / `derive`: hold or compute parameters, shrinking the search vector.
    """
    bounds = bounds or DEFAULT_BOUNDS
    fixed = fixed or {}
    derive = derive or {}
    free = _free_params(fixed, derive)
    log_bounds = [tuple(np.log(bounds[p])) for p in free]
    obj = make_objective(targets, weights=weights, sim_kwargs=sim_kwargs,
                         free=free, fixed=fixed, derive=derive)

    init = "latinhypercube"
    if x0 is not None:
        rng = np.random.default_rng(seed)
        pop = np.array([
            [rng.uniform(lo, hi) for lo, hi in log_bounds]
            for _ in range(popsize * len(free))
        ])
        pop[0] = [np.log(x0[p]) for p in free]
        init = pop

    res = differential_evolution(
        obj, log_bounds, maxiter=maxiter, popsize=popsize, seed=seed,
        init=init, tol=1e-6, mutation=(0.4, 1.0), recombination=0.7,
        polish=False, disp=disp, workers=workers,
        updating="deferred" if workers != 1 else "immediate",
    )

    x_best = res.x
    if polish_nm:
        nm = minimize(obj, x_best, method="Nelder-Mead", bounds=log_bounds,
                      options={"xatol": 1e-4, "fatol": 1e-8, "maxiter": 400})
        if nm.fun < res.fun:
            x_best = nm.x

    params_hat = _build_params(x_best, free, fixed, derive)
    return params_hat, res


def save_params(params: dict, path: str) -> None:
    with open(path, "w") as f:
        json.dump(params, f, indent=2)


def load_params(path: str) -> dict:
    with open(path) as f:
        return json.load(f)
