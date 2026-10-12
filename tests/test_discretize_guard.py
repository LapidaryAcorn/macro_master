"""
Regression test for the discretize_process() ergodic-variance guard
(DECISIONS.md, 2026-10-12): a near-random-walk component has no sensible
ergodic distribution to discretize, and silently building a chain from one
produces a grid-truncation artefact (see the guard's docstring for the
concrete incident this was written after).
"""
from __future__ import annotations

import pytest

from kmv_earnings.discretize import discretize_process
from kmv_earnings.simulate import KMV_TABLE3_PARAMS


def test_near_random_walk_component_raises():
    params = dict(KMV_TABLE3_PARAMS)
    params["beta2"] = 1e-6  # mirrors a nested-specs V3 (beta2=0) estimate
    with pytest.raises(ValueError, match="no sensible ergodic distribution"):
        discretize_process(params)


def test_infinite_variance_raises():
    params = dict(KMV_TABLE3_PARAMS)
    params["beta2"] = 0.0  # lambda*sigma^2/(2*0) = inf
    with pytest.raises(ValueError, match="no sensible ergodic distribution"):
        discretize_process(params)


def test_kmv_table3_params_still_discretize_fine():
    # sanity: the guard must not reject a normal, well-behaved process
    disc = discretize_process(KMV_TABLE3_PARAMS)
    assert "combined" in disc


def test_canonical_estimates_still_discretize_fine():
    # the pinned GER/USA canonical estimates sit right at the project's own
    # chosen ergodic-variance anchor (~0.91) and must keep working
    import json
    for path in ("output/ger/params.json", "output/usa/params.json",
                "output/fra/params.json", "output/pol/params.json"):
        params = json.load(open(path))
        disc = discretize_process(params)
        assert "combined" in disc
