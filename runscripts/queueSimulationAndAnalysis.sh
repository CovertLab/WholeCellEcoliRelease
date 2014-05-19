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

FIRST_SINGLE_ANALYSIS_JOB=""
LAST_SINGLE_ANALYSIS_JOB=""

for (( i=1; i<=$NSIMS; i++ )); do
	THIS_SIMULATION_JOB=$(qsub -v SUBMISSION_TIME=${SUBMISSION_TIME},ARRAY_ID=${i} ./runscripts/runSimulationJob.sh)
	echo THIS_SIMULATION_JOB $THIS_SIMULATION_JOB

	THIS_SINGLE_ANALYSIS_JOB=$(qsub -W depend="afterok:${THIS_SIMULATION_JOB}" -v SUBMISSION_TIME=${SUBMISSION_TIME},ARRAY_ID=${i} ./runscripts/runAnalysisSingleJob.sh)
	echo THIS_SINGLE_ANALYSIS_JOB $THIS_SINGLE_ANALYSIS_JOB

	if [ "$i" -eq 1 ]; then
		FIRST_SINGLE_ANALYSIS_JOB=${THIS_SINGLE_ANALYSIS_JOB}
	fi
	if [ "$i" -eq "$NSIMS" ]; then
		LAST_SINGLE_ANALYSIS_JOB=${THIS_SINGLE_ANALYSIS_JOB}
	fi

done

COHORT_ANALYSIS_JOB=$(qsub -W depend="afterok:${FIRST_SINGLE_ANALYSIS_JOB}:${LAST_SINGLE_ANALYSIS_JOB}" -v SUBMISSION_TIME=${SUBMISSION_TIME} ./runscripts/runAnalysisCohortJob.sh)
echo COHORT_ANALYSIS_JOB $COHORT_ANALYSIS_JOB