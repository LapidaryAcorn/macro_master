#!/usr/bin/env bash
# DECISIONS Sec 33.9, P2: 6 starts of the exact-aggregation variant on 3 lanes
# x 1 process; one committer loop commits and pushes each finished start.
set -u
cd "$(dirname "$0")/.."
L=logs/pex.log
mkdir -p logs
echo "start pex $(date)" >> $L
RD=output/pol_v4_nob/pl_exact/runs

commit_new() {
    for f in $RD/de*.json; do
        [ -e "$f" ] || continue
        if ! git ls-files --error-unmatch "$f" > /dev/null 2>&1; then
            n=$(basename "$f" .json)
            git add "$f" && git commit -q -m "Sec 33 P2: pex start $n finished

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" -- "$f" && echo "committed $n $(date)" >> $L
            git push -q >> $L 2>&1 || echo "push failed $(date)" >> $L
        fi
    done
}

python -u scripts/pol_v4.py start --variant pex --ids 1,4 --workers 1 > logs/pex_lane1.log 2>&1 &
P1=$!
python -u scripts/pol_v4.py start --variant pex --ids 2,5 --workers 1 > logs/pex_lane2.log 2>&1 &
P2=$!
python -u scripts/pol_v4.py start --variant pex --ids 3,6 --workers 1 > logs/pex_lane3.log 2>&1 &
P3=$!
while kill -0 $P1 2>/dev/null || kill -0 $P2 2>/dev/null || kill -0 $P3 2>/dev/null; do
    commit_new
    sleep 120
done
commit_new
echo "starts done $(date)" >> $L
python -u scripts/pol_v4.py finalize --variant pex > logs/pex_finalize.log 2>&1
echo "DONE_pex rc=$? $(date)" >> $L
