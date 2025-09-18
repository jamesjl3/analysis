#!/bin/bash
set -e
set -o pipefail

# env (match the release you compiled against)
export MYINSTALL=/sphenix/user/jamesj3j3/sPHENIX/install
export LD_LIBRARY_PATH="${MYINSTALL}/lib:${LD_LIBRARY_PATH:-}"
set +u
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.416
source /opt/sphenix/core/bin/setup_local.sh "${MYINSTALL}"
set -u

# lists + job index
simList=.../MatchingSubjets/dst_calo_cluster.list
truthList=.../MatchingSubjets/dst_truth_jet.list
globalList=.../MatchingSubjets/dst_global.list
g4List=.../MatchingSubjets/g4hits.list

JOB_INDEX=${1:? "Usage: $0 <job_index> [events_per_job]"}
EVENTS_PER_JOB=${2:-1000}
START_EVENT=$(( JOB_INDEX * EVENTS_PER_JOB ))
END_EVENT=$(( START_EVENT + EVENTS_PER_JOB ))

outdir=/sphenix/tg/tg01/jets/jamesj3j3/JetMatchingSubjets
mkdir -p "$outdir"
outfile="${outdir}/JetMatchingJets-w_JES-$(printf "%04d" "$JOB_INDEX").root"

# go to macro dir so relative includes work
MACRO_DIR=/sphenix/user/jamesj3j3/analysis/JS-Jet/JetValidation/macro/MatchingSubjets
cd "$MACRO_DIR"

# **single, simple call**
root.exe -l -b -q 'Fun4All_JetMatchingSubjets.C("'"$truthList"'", "'"$simList"'", "'"$globalList"'", "'"$g4List"'", "'"$outfile"'", '"$START_EVENT"', '"$END_EVENT"', '"$JOB_INDEX"')'
