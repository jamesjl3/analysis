#!/bin/bash

# Setup environment
export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/user/jamesj3j3/analysis/JS-Jet/JetValidation/macro/
export LD_LIBRARY_PATH=/sphenix/user/jamesj3j3/sPHENIX/install/lib:$LD_LIBRARY_PATH
export MYINSTALL=/sphenix/user/jamesj3j3/sPHENIX/install/

source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.502
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

printenv > env_output.txt

# Define parameters
submit=${1:-'test'}
nFile=${2:-'1000'}
nevts=${3:-'1000'}  # Events per job
minpt=${4:-'1.0'}

# File lists
simFileList=/sphenix/user/jamesj3j3/analysis/JS-Jet/JetValidation/macro/MatchingSubjets/dst_calo_cluster.list
truthList=/sphenix/user/jamesj3j3/analysis/JS-Jet/JetValidation/macro/MatchingSubjets/dst_truth_jet.list
globalList=/sphenix/user/jamesj3j3/analysis/JS-Jet/JetValidation/macro/MatchingSubjets/dst_global.list
g4hitList=/sphenix/user/jamesj3j3/analysis/JS-Jet/JetValidation/macro/MatchingSubjets/g4hits.list

# Determine number of files
if [[ $nFile -eq 0 ]]; then 
    nFile=$(wc -l < $truthList)
fi 

# Loop over files to create Condor jobs
for i in $(seq 0 $((nFile-1))); do 
    j=$(( i + 1 ))
    fname="/sphenix/user/jamesj3j3/analysis/JS-Jet/JetValidation/macro/MatchingSubjets/condor_files/condor_segment_${i}.job"
    truthf=$(sed "${j}q;d" $truthList)
    simf=$(sed "${j}q;d" $simFileList)
    globalf=$(sed "${j}q;d" $globalList)
    g4hitf=$(sed "${j}q;d" $g4hitList)

    outfile=/sphenix/tg/tg01/jets/jamesj3j3/JetMatchingSubjets/JetMatching-wJES-$(printf "%04d" $i).root

    echo "Universe  = vanilla" > $fname
    echo "Executable = /sphenix/user/jamesj3j3/analysis/JS-Jet/JetValidation/macro/MatchingSubjets/JetMatSub_Script.sh" >> $fname
    echo "Arguments = ${i} ${outfile}" >> $fname 

    echo "Output  = /sphenix/user/jamesj3j3/analysis/JS-Jet/JetValidation/macro/MatchingSubjets/condor_files/condor_${i}.out" >> $fname
    echo "Error   = /sphenix/user/jamesj3j3/analysis/JS-Jet/JetValidation/macro/MatchingSubjets/condor_files/condor_${i}.err" >> $fname
    echo "Log     = /sphenix/user/jamesj3j3/analysis/JS-Jet/JetValidation/macro/MatchingSubjets/condor_files/condor_${i}.log" >> $fname
    echo "request_memory = 128GB" >> $fname
    echo "Queue 1" >> $fname

    if [[ $submit == "submit" ]]; then 
        condor_submit $fname
    fi 

done 

echo "All jobs submitted."
