"""
TikTak global optimisation (Arnoud, Guvenen, Kleineberg) for the SMM
objective - the multistart algorithm used in the income-process estimation
literature (e.g. Guvenen et al.).

Algorithm:
  1. Global phase: evaluate the objective at n_sobol quasi-random Sobol
     points spanning the (log-)parameter box.
  2. Keep the n_local best points, sorted best-first.
  3. Local phase: run a sequence of Nelder-Mead searches. The i-th search
     starts from a convex combination of the i-th Sobol candidate and the
     best point found so far:
         start_i = (1 - w_i) * candidate_i + w_i * best_so_far,
     with w_i increasing toward 1 - early searches explore, later searches
     concentrate around the incumbent.

Uses the same objective (make_objective), bounds and common-random-numbers
simulator as estimate.py, so results are directly comparable with
differential evolution.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize
from scipy.stats import qmc

from .simulate import PARAM_ORDER
from .estimate import DEFAULT_BOUNDS, make_objective, _vec_to_params


def tiktak(
    targets: dict,
    bounds: dict | None = None,
    weights: dict | None = None,
    sim_kwargs: dict | None = None,
    n_sobol: int = 256,
    n_local: int = 10,
    local_maxfev: int = 120,
    seed: int = 11,
    disp: bool = True,
):
    """Returns (params_hat, info) with info logging the whole search."""
    bounds = bounds or DEFAULT_BOUNDS
    lb = np.array([np.log(bounds[p][0]) for p in PARAM_ORDER])
    ub = np.array([np.log(bounds[p][1]) for p in PARAM_ORDER])
    obj = make_objective(targets, weights=weights, sim_kwargs=sim_kwargs)

    # ---- global (Sobol) phase ----
    sampler = qmc.Sobol(d=len(PARAM_ORDER), scramble=True, seed=seed)
    pts = lb + qmc.Sobol.random(sampler, n_sobol) * (ub - lb)
    fvals = np.empty(n_sobol)
    for i, x in enumerate(pts):
        fvals[i] = obj(x)
        if disp and (i + 1) % 50 == 0:
            print(f"  sobol {i+1}/{n_sobol}, best so far f = {fvals[:i+1].min():.5f}")

    order = np.argsort(fvals)[:n_local]
    candidates = pts[order]
    x_best, f_best = candidates[0].copy(), float(fvals[order[0]])

    # ---- local (sequential Nelder-Mead) phase ----
    history = []
    for i, cand in enumerate(candidates):
        w = min(0.95, np.sqrt((i + 1) / (n_local + 1)))
        x0 = (1.0 - w) * cand + w * x_best
        res = minimize(obj, x0, method="Nelder-Mead",
                       options={"maxfev": local_maxfev, "xatol": 1e-4, "fatol": 1e-8})
        history.append({"start_f": float(obj(x0)), "end_f": float(res.fun)})
        if res.fun < f_best:
            x_best, f_best = res.x.copy(), float(res.fun)
        if disp:
            print(f"  local {i+1}/{n_local}: w={w:.2f}, f = {res.fun:.5f} "
                  f"(best {f_best:.5f})")

    # final polish of the incumbent
    res = minimize(obj, x_best, method="Nelder-Mead",
                   options={"maxfev": 2 * local_maxfev, "xatol": 1e-5, "fatol": 1e-9})
    if res.fun < f_best:
        x_best, f_best = res.x.copy(), float(res.fun)

    info = {"f_best": f_best, "n_evals_global": n_sobol, "history": history}
    return _vec_to_params(x_best), info
