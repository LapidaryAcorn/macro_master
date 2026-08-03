"""
Discretization of the estimated jump-drift process into a finite-state
Markov chain (the "Model Discretized" column of KMV Table 2), and export of
the chain for downstream models (KMV HANK codes / Auclert et al.
sequence-space toolkit in step 2).

Scheme (per component j):
  - bin edges power-spaced on [-zmax_mult*sigma_j, +zmax_mult*sigma_j]
    (outermost edges extended to +-inf); grid POINTS are the conditional
    means of N(0, sigma_j^2) within each bin - this keeps the size of
    discrete moves consistent with the average jump landing in the cell,
  - drift -beta_j*z: upwind finite differences (mass flows one point toward
    zero at the matching rate),
  - ADDITIVE jumps at rate lambda_j: landing probabilities integrate
    N(current state, sigma_j^2) over the bins.
Generator Q_j -> quarterly transition matrix expm(Q_j). Combined chain =
Kronecker product over the two independent components (3 x 11 = 33 states,
matching KMV: 3 transitory x 11 persistent points).

NOTE: this is a self-contained approximation. For an exact replication of
KMV's discretized column, swap in the routine from the KMV replication
package (their appendix D.1); the rest of this pipeline is unaffected.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import expm
from scipy.stats import norm

from .simulate import moments_from_panel


def component_grid_and_edges(n: int, sigma: float, zmax_mult: float = 3.0,
                             curv: float = 1.8):
    x = np.linspace(-1.0, 1.0, n + 1)
    edges = np.sign(x) * np.abs(x) ** curv * zmax_mult * sigma
    edges[0], edges[-1] = -np.inf, np.inf
    pts = np.empty(n)
    for k in range(n):
        a, b = edges[k], edges[k + 1]
        mass = norm.cdf(b, scale=sigma) - norm.cdf(a, scale=sigma)
        pa = norm.pdf(a / sigma) if np.isfinite(a) else 0.0
        pb = norm.pdf(b / sigma) if np.isfinite(b) else 0.0
        pts[k] = sigma * (pa - pb) / max(mass, 1e-12)
    return pts, edges


def component_generator(grid: np.ndarray, edges: np.ndarray, lam: float,
                        beta: float, sigma: float) -> np.ndarray:
    n = len(grid)
    Q = np.zeros((n, n))
    for k in range(n):                      # drift toward zero (upwind)
        g = grid[k]
        if g > 0 and k > 0:
            r = beta * g / (grid[k] - grid[k - 1])
            Q[k, k - 1] += r
            Q[k, k] -= r
        elif g < 0 and k < n - 1:
            r = beta * (-g) / (grid[k + 1] - grid[k])
            Q[k, k + 1] += r
            Q[k, k] -= r
    for k in range(n):                      # additive jumps
        w = np.diff(norm.cdf(edges, loc=grid[k], scale=sigma))
        Q[k, :] += lam * w
        Q[k, k] -= lam
    return Q


def stationary_dist(P: np.ndarray) -> np.ndarray:
    vals, vecs = np.linalg.eig(P.T)
    k = int(np.argmin(np.abs(vals - 1.0)))
    pi = np.abs(np.real(vecs[:, k]))
    return pi / pi.sum()


def discretize_process(params: dict, n1: int = 3, n2: int = 11,
                       zmax_mult: float = 3.0, curv: float = 1.8):
    """n1 = transitory grid points, n2 = persistent (3 x 11 = 33 as in KMV)."""
    out = {}
    comps = [
        ("1", params["lambda1"], params["beta1"], params["sigma1"], n1),
        ("2", params["lambda2"], params["beta2"], params["sigma2"], n2),
    ]
    for name, lam, beta, sig, n in comps:
        grid, edges = component_grid_and_edges(n, sig, zmax_mult, curv)
        Q = component_generator(grid, edges, lam, beta, sig)
        P = expm(Q)  # quarterly transition matrix
        out[name] = {"grid": grid, "Q": Q, "P": P, "pi": stationary_dist(P)}

    g1, g2 = out["1"]["grid"], out["2"]["grid"]
    zgrid = (g1[:, None] + g2[None, :]).ravel()
    P = np.kron(out["1"]["P"], out["2"]["P"])
    pi = np.kron(out["1"]["pi"], out["2"]["pi"])
    e = np.exp(zgrid)
    e = e / (pi @ e)  # normalise mean earnings to 1
    out["combined"] = {"zgrid": zgrid, "egrid": e, "P": P, "pi": pi}
    return out


def simulate_discrete_panel(disc: dict, n_workers: int = 50_000,
                            n_years_keep: int = 36, lifecycle: bool = True,
                            seed: int = 4321) -> np.ndarray:
    """
    Simulate the chain at its native quarterly frequency (the frequency at
    which the HANK model uses it), aggregate the flow exp(z) to annual.
    lifecycle=True starts all workers at the grid point closest to zero,
    mirroring the continuous lifecycle design.
    """
    rng = np.random.default_rng(seed)
    logy = np.empty((n_workers, n_years_keep))
    states, cums = {}, {}
    for name in ("1", "2"):
        grid, pi = disc[name]["grid"], disc[name]["pi"]
        if lifecycle:
            states[name] = np.full(n_workers, int(np.argmin(np.abs(grid))))
        else:
            states[name] = rng.choice(len(pi), size=n_workers, p=pi)
        cums[name] = np.cumsum(disc[name]["P"], axis=1)

    annual_flow = np.zeros(n_workers)
    year_idx = 0
    for q in range(4 * n_years_keep):
        for name in ("1", "2"):
            u = rng.random(n_workers)
            states[name] = (u[:, None] > cums[name][states[name]]).sum(axis=1)
        z = disc["1"]["grid"][states["1"]] + disc["2"]["grid"][states["2"]]
        annual_flow += np.exp(z)
        if q % 4 == 3:
            logy[:, year_idx] = np.log(annual_flow)
            annual_flow[:] = 0.0
            year_idx += 1
    return logy


def discretized_moments(params: dict, disc: dict | None = None, **sim_kwargs) -> dict:
    disc = disc or discretize_process(params)
    return moments_from_panel(simulate_discrete_panel(disc, **sim_kwargs))


def export_chain(disc: dict, outdir: str, prefix: str = "income_process") -> None:
    """
    Write the combined chain to text files for step 2:
      *_grid.txt (earnings levels, mean 1), *_zgrid.txt (log grid),
      *_P.txt (quarterly 33x33 transition matrix), *_pi.txt (stationary dist).
    """
    import os
    c = disc["combined"]
    np.savetxt(os.path.join(outdir, f"{prefix}_grid.txt"), c["egrid"])
    np.savetxt(os.path.join(outdir, f"{prefix}_zgrid.txt"), c["zgrid"])
    np.savetxt(os.path.join(outdir, f"{prefix}_P.txt"), c["P"])
    np.savetxt(os.path.join(outdir, f"{prefix}_pi.txt"), c["pi"])
