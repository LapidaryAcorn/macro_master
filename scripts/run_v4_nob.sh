#!/usr/bin/env bash
# Sec 33 (v4_nob): 12 starts for one side, 3 lanes x 1 process (max 3 processes:
# the workstation powers off under sustained full load, Sec 32.9), then finalize.
# Usage: bash scripts/run_v4_nob.sh nob|kmv_nob
# Restartable: finished starts (runs/deNN.json) are skipped.
set -u
cd "$(dirname "$0")/.."
V=$1
L=logs/v4_nob_$V.log
mkdir -p logs
echo "start $V $(date)" >> $L
python -u scripts/pol_v4.py start --variant $V --ids 1,4,7,10 --workers 1 > logs/v4_nob_${V}_lane1.log 2>&1 &
python -u scripts/pol_v4.py start --variant $V --ids 2,5,8,11 --workers 1 > logs/v4_nob_${V}_lane2.log 2>&1 &
python -u scripts/pol_v4.py start --variant $V --ids 3,6,9,12 --workers 1 > logs/v4_nob_${V}_lane3.log 2>&1 &
wait
echo "starts done $(date)" >> $L
python -u scripts/pol_v4.py finalize --variant $V > logs/v4_nob_${V}_finalize.log 2>&1
echo "DONE_$V rc=$? $(date)" >> $L
