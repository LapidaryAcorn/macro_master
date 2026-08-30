"""
Add KMV-style stochastic death / perpetual youth to the ARS one-asset HA block.

KMV (2018) footnote 15: households die at rate zeta = 1/180 per quarter (45-year
expected working life) and are replaced by newborns with ZERO wealth; this is the
mechanism that puts enough households at zero liquid wealth to generate the
aggregate MPC. ARS's one-asset code is infinitely-lived and has no such device
(confirmed: no death/mortality/rebirth anywhere in annual-review/).

Two modifications, nothing else:
  1. Backward (discount): the household values continuation only if alive, so the
     effective discount is beta * (1 - zeta).  No bequest motive.
  2. Forward (distribution): each period a fraction zeta of the mass at every
     (beta, e, a) is moved to asset index 0 (a = 0 = a_grid[0]).  Newborns keep
     their (beta, e) state -- but since the cross-sectional (beta, e) distribution
     IS the stationary one at all times (both are exogenous with stationary
     pi_b, pi_e), newborns are drawn from the ergodic (beta, e) distribution in
     the cross-section, exactly as in KMV.  So earnings dispersion is unchanged:
     the chain's ergodic var(log e) still governs the cross-section.

Implemented by subclassing SSJ's HetBlock and overriding make_endog_law_of_motion
with a death-modified policy lottery, so it flows through both steady_state and
jacobian.
"""

from __future__ import annotations
import numpy as np
import sequence_jacobian as sj
from sequence_jacobian.blocks.het_block import HetBlock
from sequence_jacobian.blocks.support.het_support import (
    PolicyLottery1D, ForwardShockablePolicyLottery1D, lottery_1d)


# ---- death-modified policy lottery ----------------------------------------
class DeathLottery1D(PolicyLottery1D):
    def __init__(self, i, pi, grid, zeta):
        super().__init__(i, pi, grid)
        self.zeta = float(zeta)

    def forward(self, D):
        Dsurv = super().forward(D)                      # normal asset lottery
        out = (1.0 - self.zeta) * Dsurv.reshape(self.flatshape)
        out[:, 0] += self.zeta * D.reshape(self.flatshape).sum(axis=1)   # newborns at a=0
        return out.reshape(self.shape)

    def expectation(self, X):
        Xsurv = super().expectation(X).reshape(self.flatshape)
        Xflat = X.reshape(self.flatshape)
        newborn = np.repeat(Xflat[:, :1], Xflat.shape[1], axis=1)        # value at a=0, same (beta,e)
        return ((1.0 - self.zeta) * Xsurv + self.zeta * newborn).reshape(self.shape)

    def forward_shockable(self, Dss):
        return ForwardShockableDeathLottery1D(
            self.i.reshape(self.shape), self.pi.reshape(self.shape), self.grid, Dss, self.zeta)


class ForwardShockableDeathLottery1D(ForwardShockablePolicyLottery1D):
    def __init__(self, i, pi, grid, Dss, zeta):
        super().__init__(i, pi, grid, Dss)
        self.zeta = float(zeta)

    def forward_shock(self, da):
        # newborn mass zeta * (row total) does not depend on the policy shock da,
        # so d(forward)/d(policy) is just (1 - zeta) times the normal response
        return (1.0 - self.zeta) * super().forward_shock(da)


# ---- HetBlock that uses it ----------------------------------------------
class HetBlockDeath(HetBlock):
    zeta_death = 0.0                      # set on the instance before use

    def make_endog_law_of_motion(self, d, monotonic=False):
        lot = lottery_1d(d[self.policy[0]], d[self.policy[0] + '_grid'], monotonic)
        return DeathLottery1D(lot.i.reshape(lot.shape), lot.pi.reshape(lot.shape),
                              lot.grid, self.zeta_death)


def _hh_raw_death(Va_p, a_grid, y, r, beta, eis, zeta_d, annuity):
    """ARS hh_raw with Blanchard-Yaari perpetual youth.
      - continuation discounted by survival: beta -> beta * (1 - zeta_d)
      - annuity=1: survivors earn the deceased's assets, gross return
        (1+r)/(1-zeta_d); aggregate wealth conserved (standard B-Y).
        annuity=0: no annuities, gross return (1+r); the deceased's wealth
        leaks out of the household sector (must be rebated elsewhere in GE).
    """
    Rg = (1.0 + r) / (1.0 - zeta_d * annuity)
    uc_nextgrid = (beta[:, np.newaxis] * (1.0 - zeta_d)) * Va_p
    c_nextgrid = uc_nextgrid ** (-eis)
    coh = Rg * a_grid[np.newaxis, :] + y[:, np.newaxis]
    a = sj.interpolate.interpolate_y(c_nextgrid + a_grid, coh, a_grid)
    sj.misc.setmin(a, a_grid[0])
    c = coh - a
    Va = Rg * c ** (-1 / eis)
    return Va, a, c


def make_hh_death(zeta_death=1/180):
    hb = HetBlockDeath(_hh_raw_death, exogenous='Pi', policy='a', backward='Va',
                       backward_init=sj.hetblocks.hh_sim.hh_init)
    hb.zeta_death = float(zeta_death)
    return hb
