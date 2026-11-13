#!/usr/bin/env bash
# Sec 32, phase 1: 21 DE starts of the v4 baseline, 3 lanes x 4 DE workers
# (restartable: finished starts are skipped).
set -u
cd "$(dirname "$0")/.."
S=scripts/pol_v4.py
mkdir -p logs
echo "phase1 start $(date)" >> logs/v4_phase1.log
python -u $S start --variant base --ids 1,4,7,10,13,16,19 --workers 4 >> logs/v4_lane1.log 2>&1 &
python -u $S start --variant base --ids 2,5,8,11,14,17,20 --workers 4 >> logs/v4_lane2.log 2>&1 &
python -u $S start --variant base --ids 3,6,9,12,15,18,21 --workers 4 >> logs/v4_lane3.log 2>&1 &
wait
echo "PHASE1_DONE $(date)" >> logs/v4_phase1.log
