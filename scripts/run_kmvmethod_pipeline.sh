#!/usr/bin/env bash
# Sec 31.3-31.4 driver (restartable: finished starts are skipped).
#  A: 3 kmv DE seeds (4 workers each)  ||  B: small-budget decomposition DE (4 workers)
#  after A: kmv TikTak (6 workers) + kmv DFLS starts; then finalize everything.
set -u
cd "$(dirname "$0")/.."
S=scripts/pol_v3_kmvmethod.py
EL="v3 reset euler16 ergodic gap4 kurtw box nsim5000"
echo "start $(date)"
( for v in $EL; do python -u $S --variant $v --run de --seeds 7 --small --workers 4; done ) >> logs/kmvm_decomp_de.log 2>&1 &
PB=$!
for s in 7 101 202; do python -u $S --variant kmv --run de --seeds $s --workers 4 > logs/kmvm_pol_de$s.log 2>&1 & done
for j in $(jobs -p); do [ "$j" != "$PB" ] && wait $j; done
echo "kmv DE done $(date)"
python -u $S --variant kmv --run dfls > logs/kmvm_pol_dfls.log 2>&1 &
python -u $S --variant kmv --run dfls_v3 >> logs/kmvm_pol_dfls2.log 2>&1 &
python -u $S --variant kmv --run tiktak --workers 6 > logs/kmvm_pol_tiktak.log 2>&1
wait
echo "finalize $(date)"
python -u $S --variant kmv --finalize > logs/kmvm_finalize.log 2>&1
for v in $EL dfls; do python -u $S --variant $v --finalize >> logs/kmvm_finalize.log 2>&1; done
echo "PIPELINE_DONE $(date)"
