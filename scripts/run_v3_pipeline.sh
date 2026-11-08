#!/usr/bin/env bash
# Sec 26-28: run the Poland v3 steps that depend on the canonical v3 estimate,
# one after another, once all six v3 starts (5 DE seeds + TikTak) are on disk.
# Each step writes its own log under output/pol_v3/logs/. Restartable: every
# script skips work whose output already exists.
set -u
cd "$(dirname "$0")/.."
L=output/pol_v3/logs
step() { echo "[$(date +%H:%M)] START $1" >> $L/pipeline.log; shift_name=$1; shift
         "$@" > $L/$shift_name.log 2>&1; echo "[$(date +%H:%M)] END $shift_name rc=$?" >> $L/pipeline.log; }
until [ "$(ls output/pol_v3/runs/de_seed*.json output/pol_v3/runs/tiktak.json 2>/dev/null | wc -l)" -ge 6 ]; do sleep 60; done
step finalize   python -u scripts/pol_v3_estimate.py --version v3 --finalize
step smm_se     python -u scripts/pol_v3_smm_se.py --workers 10
step figures    python -u scripts/pol_v3_figures.py
step vs_kmv     python -u scripts/pol_v3_vs_kmv.py table
step bootstrap  python -u scripts/pol_v3_bootstrap.py --workers 10
step profile    python -u scripts/pol_v3_robustness.py profile --workers 10
step recovery   python -u scripts/pol_v3_robustness.py recovery --workers 10
step nested     python -u scripts/pol_v3_robustness.py nested --workers 10
step robust_sum python -u scripts/pol_v3_robustness.py summary
step v3_end_de  python -u scripts/pol_v3_estimate.py --version v3_end --only de --workers 10
until [ -f output/pol_v3_end/runs/tiktak.json ]; do sleep 60; done
step v3_end_fin python -u scripts/pol_v3_estimate.py --version v3_end --finalize
step grid_size  python -u scripts/pol_v3_vs_kmv.py grid
echo "[$(date +%H:%M)] PIPELINE DONE" >> $L/pipeline.log
