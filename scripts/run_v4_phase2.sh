#!/usr/bin/env bash
# Sec 32, phase 2 (after phase 1): finalize, SEs, bootstrap, profile, recovery,
# robustness variants, figures, HANK steady state, nested specs.
# Load capped at 6 processes after a hard power-off under 12 (2026-11-13).
# Restartable: every item skips what already exists.
set -u
cd "$(dirname "$0")/.."
S=scripts/pol_v4.py; U=scripts/pol_v4_uncertainty.py
L=logs/v4_phase2.log
until grep -q PHASE1_DONE logs/v4_phase1.log 2>/dev/null; do sleep 120; done
echo "phase2 start $(date)" >> $L
[ -f output/pol_v4/estimate/report.json ] || python -u $S finalize --variant base > logs/v4_finalize_base.log 2>&1
echo "finalize done $(date)" >> $L
[ -f output/pol_v4/smm_se/smm_se.json ] || python -u $U smm_se --workers 6 > logs/v4_smm_se.log 2>&1
echo "smm_se done $(date)" >> $L
python -u $U boot --workers 6 > logs/v4_boot.log 2>&1
echo "boot done $(date)" >> $L
python -u $U profile --workers 6 --lanes 2 > logs/v4_profile.log 2>&1
echo "profile done $(date)" >> $L
python -u $U recovery --workers 6 --lanes 2 > logs/v4_recovery.log 2>&1
echo "recovery done $(date)" >> $L
for v in full5 nob1 eqw; do python -u $S start --variant $v --ids 1,2,3 --workers 2 > logs/v4_robust_$v.log 2>&1 & done
wait
for v in full5 nob1 eqw; do python -u $S finalize --variant $v >> logs/v4_robust_$v.log 2>&1; done
echo "robust done $(date)" >> $L
python -u scripts/pol_v4_figures.py > logs/v4_figures.log 2>&1
python -u hank_step2/v4_steady_state.py > logs/v4_hank.log 2>&1
echo "figures+hank done $(date)" >> $L
python -u $U nested --workers 6 --lanes 2 > logs/v4_nested.log 2>&1
python -u $U summary > logs/v4_summary.log 2>&1
echo "PHASE2_DONE $(date)" >> $L
