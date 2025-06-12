#!/bin/bash

# Setup environment
export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/user/jamesj3j3/analysis/JS-Jet/JetValidation/macro/
export LD_LIBRARY_PATH=/sphenix/user/jamesj3j3/sPHENIX/install/lib:$LD_LIBRARY_PATH
export MYINSTALL=/sphenix/user/jamesj3j3/sPHENIX/install/

source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.416
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

printenv > env_output.txt

# Define total events and events per job
TOTAL_EVENTS=10000000
EVENTS_PER_JOB=10000  # Adjust as needed

# Get job index from Condor
JOB_INDEX=$1
outfile=$2
echo "JOB_INDEX: $JOB_INDEX"

# Calculate start and end event number for this job
START_EVENT=$((JOB_INDEX * EVENTS_PER_JOB))
END_EVENT=$((START_EVENT + EVENTS_PER_JOB))

outfile=/sphenix/tg/tg01/jets/jamesj3j3/JetMatchingSubjets/JetMatchingSubjets-$(printf "%04d" $JOB_INDEX).root
#outfile=/sphenix/tg/tg01/jets/jamesj3j3/JetMatchingSubjets/JetMatchingSubjets-$(printf "%04d" $JOB_INDEX).root
echo "Generated outfile: $outfile"

# File lists
simFileList=/sphenix/user/jamesj3j3/analysis/JS-Jet/JetValidation/macro/JetCorrection/dst_calo_cluster.list
truthList=/sphenix/user/jamesj3j3/analysis/JS-Jet/JetValidation/macro/JetCorrection/dst_truth_jet.list
globalList=/sphenix/user/jamesj3j3/analysis/JS-Jet/JetValidation/macro/JetCorrection/dst_global.list
g4hitList=/sphenix/user/jamesj3j3/analysis/JS-Jet/JetValidation/macro/JetCorrection/dst_truth.list

# Get the corresponding files for this job index
truthFile=$(sed "$((JOB_INDEX + 1))q;d" $truthList)
echo "Selected truthFile: $truthFile"
simFile=$(sed "$((JOB_INDEX + 1))q;d" $simFileList)
globalFile=$(sed "$((JOB_INDEX + 1))q;d" $globalList)
g4hitFile=$(sed "$((JOB_INDEX + 1))q;d" $g4hitList)

echo "Processing events from $START_EVENT to $END_EVENT"
echo "truthFile: $truthFile"
echo "simFile: $simFile"
echo "globalFile: $globalFile"
echo "g4hitFile: $g4hitFile"
echo "outfile: $outfile"

# Run the Fun4All macro with the selected files
#valgrind
root.exe -q -b Fun4All_JetMatchingSubjets.C\(\"$truthFile\",\"$simFile\",\"$globalFile\",\"$g4hitFile\",\"$outfile\",$START_EVENT,$END_EVENT\)
echo "root.exe -q -b Fun4All_JetMatchingSubjets.C\(\"$truthFile\",\"$simFile\",\"$globalFile\",\"$g4hitFile\",\"$outfile\",$START_EVENT,$END_EVENT\)"
echo "Job $JOB_INDEX completed."
