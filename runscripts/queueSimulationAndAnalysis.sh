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
KB_FIT="${KB_DIR}/KnowledgeBase_Most_Fit.cPickle"
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

if [ -z "${VARIANT}" ]; then
	VARIANT="wildtype"
fi
if [ -z "${FIRST_VARIANT_INDEX}" ]; then
	FIRST_VARIANT_INDEX="0"
fi
if [ -z "${LAST_VARIANT_INDEX}" ]; then
	LAST_VARIANT_INDEX="0"
fi
if [ "${LAST_VARIANT_INDEX}" = "-1" ]; then
	LAST_VARIANT_INDEX="$(($(PYTHONPATH="${PWD}:${PYTHONPATH}" python2.7 ./runscripts/getNumVariants.py "${VARIANT}" "${KB_FIT}") - 1))"
fi
echo LAST_VARIANT_INDEX $LAST_VARIANT_INDEX

for (( i=${FIRST_VARIANT_INDEX}; i<=${LAST_VARIANT_INDEX}; i++ )); do
	VARIANT_DIR="out/${SUBMISSION_TIME}/${VARIANT}_$(printf "%06d" ${i})"
	VARIANT_KB_DIR="${VARIANT_DIR}/kb"
	VARIANT_METADATA_DIR="${VARIANT_DIR}/metadata"

	mkdir -p "${VARIANT_DIR}"
	mkdir -p "${VARIANT_KB_DIR}"
	mkdir -p "${VARIANT_METADATA_DIR}"

	PYTHONPATH="${PWD}:${PYTHONPATH}" python2.7 ./runscripts/createVariantKb.py "${VARIANT}" "${i}" "${KB_FIT}" "${VARIANT_KB_DIR}/KnowledgeBase_Modified.cPickle" "${VARIANT_METADATA_DIR}"


	for (( j=1; j<=$NSIMS; j++ )); do
		THIS_SIMULATION_JOB=$(qsub -v SUBMISSION_TIME=${SUBMISSION_TIME},\
VARIANT=${VARIANT},\
VARIANT_ID=${i},\
SIM_ID=${j},\
WC_LENGTHSEC=${WC_LENGTHSEC},\
WC_LOGTOSHELL=${WC_LOGTOSHELL},\
WC_LOGTODISKEVERY=${WC_LOGTODISKEVERY} ./runscripts/runSimulationJob.sh)

		THIS_SINGLE_ANALYSIS_JOB=$(qsub -W depend="afterok:${THIS_SIMULATION_JOB}" -v SUBMISSION_TIME=${SUBMISSION_TIME},VARIANT=${VARIANT},VARIANT_ID=${i},SIM_ID=${j} ./runscripts/runAnalysisSingleJob.sh)
		echo THIS_SINGLE_ANALYSIS_JOB $THIS_SINGLE_ANALYSIS_JOB

		if [ "$j" -eq 1 ]; then
			FIRST_SINGLE_ANALYSIS_JOB=${THIS_SINGLE_ANALYSIS_JOB}
		fi
		if [ "$j" -eq "$NSIMS" ]; then
			LAST_SINGLE_ANALYSIS_JOB=${THIS_SINGLE_ANALYSIS_JOB}
		fi

	done

	# COHORT_ANALYSIS_JOB=$(qsub -W depend="afterok:${FIRST_SINGLE_ANALYSIS_JOB}:${LAST_SINGLE_ANALYSIS_JOB}" -v SUBMISSION_TIME=${SUBMISSION_TIME} ./runscripts/runAnalysisCohortJob.sh)
	# echo COHORT_ANALYSIS_JOB $COHORT_ANALYSIS_JOB

	KB_COMPRESSION_JOB=$(qsub -W depend="afterany:${FIRST_SINGLE_ANALYSIS_JOB}:${LAST_SINGLE_ANALYSIS_JOB}" -v SUBMISSION_TIME=${SUBMISSION_TIME},VARIANT=${VARIANT},VARIANT_ID=${i} ./runscripts/runCompressKbsJob.sh)
	echo KB_COMPRESSION_JOB $KB_COMPRESSION_JOB

done

KB_COMPRESSION_JOB=$(qsub -v SUBMISSION_TIME=${SUBMISSION_TIME} ./runscripts/runCompressKbsJob.sh)
echo KB_COMPRESSION_JOB $KB_COMPRESSION_JOB
