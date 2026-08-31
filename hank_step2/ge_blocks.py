"""
GE (non-household) blocks for the ARS one-asset HANK, copied VERBATIM from
`Annual Review main.ipynb` section 1 (they live in the notebook, not
household.py). Imported by baseline.py and swap_chain.py so both use exactly the
same block definitions.
"""
import sequence_jacobian as sj


@sj.simple
def production(Y, mu):
    N = Y
    wN = Y / mu
    div = Y - wN
    return N, wN, div


@sj.simple
def real_ST_bonds(r_ante):
    r_post_bonds = r_ante(-1)
    qbond = 1/(1 + r_ante)
    return r_post_bonds, qbond


@sj.simple
def nom_ST_bonds(r_ante, pi):
    pi_e = pi(1)
    r_post_bonds = (1 + r_ante(-1)) * (1 + pi_e(-1))/(1 + pi) - 1
    qbond = 1/(1 + r_ante)
    return r_post_bonds, qbond


@sj.simple
def fiscal(B, r_post_bonds, G, Y, wN, div, tax_rate_shock, T_shock, T_rule_coeff, qbond):
    T = (1 + r_post_bonds) * B(-1) + G - B
    tax_rate = T / Y
    wN_aftertax = (1 - tax_rate) * wN
    div_aftertax = (1 - tax_rate) * div
    tax_rate_res = tax_rate - tax_rate.ss - tax_rate_shock
    T_res = T - T.ss - T_rule_coeff * (B(-1) - B.ss) - T_shock
    constant_owed_res = B/qbond - B.ss/qbond.ss
    return T, tax_rate, wN_aftertax, div_aftertax, tax_rate_res, T_res, constant_owed_res


@sj.solved(unknowns={'p': (0.001, 50)}, targets=['cap_cond'])
def capitalization(div_aftertax, r_ante, p):
    cap_cond = p - (div_aftertax(+1) + p(+1))/(1 + r_ante)
    r_post_equity = (div_aftertax + p)/p(-1) - 1
    return r_post_equity, cap_cond


@sj.simple
def ex_post_r(r_post_equity, r_post_bonds, B, p):
    r = (p(-1) * r_post_equity + B(-1) * r_post_bonds) / (p(-1) + B(-1))
    return r


@sj.simple
def nkpc(pi, kappa, N, Y, G, beta_ave, vscale, frisch, tax_rate, mu, eis):
    wedge = vscale*N**(1/frisch) - (1-tax_rate) / mu * (Y-G)**(-1/eis)
    nkpc_res = kappa * wedge + beta_ave * pi(+1) - pi
    return nkpc_res


@sj.simple
def taylor_rule(rstar, pi, phi):
    i = rstar + phi * pi
    r_ante = i - pi(+1)
    return r_ante


@sj.simple
def mkt_clearing(A, Y, C, p, B, G):
    asset_mkt = A - p - B
    goods_mkt = C + G - Y
    return asset_mkt, goods_mkt
