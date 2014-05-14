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

FIRST_ANALYSIS_JOB=""
LAST_ANALYSIS_JOB=""

for (( i=1; i<=$NSIMS; i++ )); do
	THIS_SIMULATION_JOB=$(qsub -v SUBMISSION_TIME=${SUBMISSION_TIME},ARRAY_ID=${i} ./runscripts/runSimulationJob.sh)
	echo THIS_SIMULATION_JOB $THIS_SIMULATION_JOB

	THIS_ANALYSIS_JOB=$(qsub -W depend="afterok:${THIS_SIMULATION_JOB}" -v SUBMISSION_TIME=${SUBMISSION_TIME},ARRAY_ID=${i} ./runscripts/runAnalysisSingleJob.sh)
	echo THIS_ANALYSIS_JOB $THIS_ANALYSIS_JOB

	if [ "$i" -eq 1 ]; then
		FIRST_ANALYSIS_JOB=${THIS_ANALYSIS_JOB}
	fi
	if [ "$i" -eq "$NSIMS" ]; then
		LAST_ANALYSIS_JOB=${THIS_ANALYSIS_JOB}
	fi

done

echo qsub -W depend="${FIRST_ANALYSIS_JOB}:${LAST_ANALYSIS_JOB}" ./runAnalysisCohortJob.sh
# ANALYSIS_COHORT_JOB=$(qsub -W depend="${FIRST_ANALYSIS_JOB}:${LAST_ANALYSIS_JOB}" ./runAnalysisCohortJob.sh)
# echo ANALYSIS_COHORT_JOB $ANALYSIS_COHORT_JOB