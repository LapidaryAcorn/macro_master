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
from .estimate import (DEFAULT_BOUNDS, make_objective, _build_params,
                       _free_params)


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
    workers: int = 1,
    fixed: dict | None = None,
    derive: dict | None = None,
    x_start: dict | None = None,
):
    """Returns (params_hat, info) with info logging the whole search.

    workers > 1 parallelises the global (Sobol) phase over a process pool; the
    local phase is inherently sequential (each start depends on the incumbent).
    `fixed` / `derive`: hold or compute parameters, shrinking the search vector
    (same semantics as estimate.estimate).
    `x_start`: optional warm-start parameter dict, added as an extra candidate
    and used as the initial incumbent if it beats the best Sobol point - so the
    DE-vs-TikTak comparison can be run with symmetric starting information.
    """
    bounds = bounds or DEFAULT_BOUNDS
    fixed = fixed or {}
    derive = derive or {}
    free = _free_params(fixed, derive)
    lb = np.array([np.log(bounds[p][0]) for p in free])
    ub = np.array([np.log(bounds[p][1]) for p in free])
    nm_bounds = list(zip(lb, ub))
    obj = make_objective(targets, weights=weights, sim_kwargs=sim_kwargs,
                         free=free, fixed=fixed, derive=derive)

    # ---- global (Sobol) phase ----
    sampler = qmc.Sobol(d=len(free), scramble=True, seed=seed)
    pts = lb + qmc.Sobol.random(sampler, n_sobol) * (ub - lb)
    if workers and workers != 1:
        import multiprocessing as mp
        with mp.Pool(None if workers < 0 else workers) as pool:
            fvals = np.array(pool.map(obj, list(pts)))
    else:
        fvals = np.empty(n_sobol)
        for i, x in enumerate(pts):
            fvals[i] = obj(x)
            if disp and (i + 1) % 50 == 0:
                print(f"  sobol {i+1}/{n_sobol}, best so far f = {fvals[:i+1].min():.5f}")

    order = np.argsort(fvals)[:n_local]
    candidates = pts[order]
    x_best, f_best = candidates[0].copy(), float(fvals[order[0]])

    if x_start is not None:
        xs = np.clip(np.array([np.log(x_start[p]) for p in free]), lb, ub)
        fs = float(obj(xs))
        if disp:
            print(f"  warm start f = {fs:.5f}  (best Sobol {f_best:.5f})")
        candidates = np.vstack([xs, candidates])
        if fs < f_best:
            x_best, f_best = xs.copy(), fs

    # ---- local (sequential Nelder-Mead) phase ----
    history = []
    for i, cand in enumerate(candidates):
        w = min(0.95, np.sqrt((i + 1) / (n_local + 1)))
        x0 = (1.0 - w) * cand + w * x_best
        res = minimize(obj, x0, method="Nelder-Mead", bounds=nm_bounds,
                       options={"maxfev": local_maxfev, "xatol": 1e-4, "fatol": 1e-8})
        history.append({"start_f": float(obj(x0)), "end_f": float(res.fun)})
        if res.fun < f_best:
            x_best, f_best = res.x.copy(), float(res.fun)
        if disp:
            print(f"  local {i+1}/{n_local}: w={w:.2f}, f = {res.fun:.5f} "
                  f"(best {f_best:.5f})")

    # final polish of the incumbent
    res = minimize(obj, x_best, method="Nelder-Mead", bounds=nm_bounds,
                   options={"maxfev": 2 * local_maxfev, "xatol": 1e-5, "fatol": 1e-9})
    if res.fun < f_best:
        x_best, f_best = res.x.copy(), float(res.fun)

    info = {"f_best": f_best, "n_evals_global": n_sobol, "history": history}
    return _build_params(x_best, free, fixed, derive), info
