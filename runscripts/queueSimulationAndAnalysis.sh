#!/bin/bash

if [ "$(basename $PWD)" != "wcEcoli" ]; then
	echo "Must execute this script from the wcEcoli 'root' directory (e.g. ~/Documents/wcEcoli)." >&2
	exit 1
fi

if [ -z "$1" ]; then
	NSIMS="1"
else
	NSIMS="$1"
fi

SUBMISSION_TIME=$(date "+%Y%m%d.%H%M%S.%N")
KB_DIR="out/${SUBMISSION_TIME}/kb"
KB_FIT="${KB_DIR}/KnowledgeBase_Fit.cPickle"
METADATA_DIR="out/${SUBMISSION_TIME}/metadata"

##### Create metadata directory #####
mkdir -p "${METADATA_DIR}"

##### Save metadata #####
echo "Adding simulation metadata"

# Git hash
git rev-parse HEAD > "${METADATA_DIR}/git_hash"

# Git branch
git symbolic-ref --short HEAD > "${METADATA_DIR}/git_branch"

# Git diff
git diff > "${METADATA_DIR}/git_diff"

# Description
echo "${DESC}" > "${METADATA_DIR}/description"

##### Create knowledgebases (unfit and fit) #####
make justKb FIXTURES_KBDIR="${KB_DIR}"

FIRST_SINGLE_ANALYSIS_JOB=""
LAST_SINGLE_ANALYSIS_JOB=""

for (( i=1; i<=$NSIMS; i++ )); do
	THIS_SIMULATION_JOB=$(qsub -v SUBMISSION_TIME=${SUBMISSION_TIME},\
ARRAY_ID=${i},\
WC_LENGTHSEC=${WC_LENGTHSEC},\
WC_LOGTOSHELL=${WC_LOGTOSHELL},\
WC_LOGTODISKEVERY=${WC_LOGTODISKEVERY} ./runscripts/runSimulationJob.sh)

	THIS_SINGLE_ANALYSIS_JOB=$(qsub -W depend="afterok:${THIS_SIMULATION_JOB}" -v SUBMISSION_TIME=${SUBMISSION_TIME},ARRAY_ID=${i} ./runscripts/runAnalysisSingleJob.sh)
	echo THIS_SINGLE_ANALYSIS_JOB $THIS_SINGLE_ANALYSIS_JOB

	if [ "$i" -eq 1 ]; then
		FIRST_SINGLE_ANALYSIS_JOB=${THIS_SINGLE_ANALYSIS_JOB}
	fi
	if [ "$i" -eq "$NSIMS" ]; then
		LAST_SINGLE_ANALYSIS_JOB=${THIS_SINGLE_ANALYSIS_JOB}
	fi

done

# COHORT_ANALYSIS_JOB=$(qsub -W depend="afterok:${FIRST_SINGLE_ANALYSIS_JOB}:${LAST_SINGLE_ANALYSIS_JOB}" -v SUBMISSION_TIME=${SUBMISSION_TIME} ./runscripts/runAnalysisCohortJob.sh)
# echo COHORT_ANALYSIS_JOB $COHORT_ANALYSIS_JOB

KB_COMPRESSION_JOB=$(qsub -W depend="afterany:${FIRST_SINGLE_ANALYSIS_JOB}:${LAST_SINGLE_ANALYSIS_JOB}" -v SUBMISSION_TIME=${SUBMISSION_TIME} ./runscripts/runCompressKbsJob.sh)
echo KB_COMPRESSION_JOB $KB_COMPRESSION_JOB