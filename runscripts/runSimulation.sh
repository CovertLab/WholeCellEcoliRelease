#!/bin/bash

trap "exit 130" 2 9 15

SUBMISSION_TIME=$(date "+%Y%m%d.%H%M%S.%N")
WC_SEED=${WC_SEED:-0}
SEED_DIR=$(printf "%06d" $(($WC_SEED)))
SIM_OUT_DATA_DIR="out/${SUBMISSION_TIME}/wildtype_000000/${SEED_DIR}/simOut"
KB_DIR="out/${SUBMISSION_TIME}/kb"
KB_FIT="${KB_DIR}/KnowledgeBase_Most_Fit.cPickle"
METADATA_DIR="out/${SUBMISSION_TIME}/metadata"
VARIANT_KB_DIR="out/${SUBMISSION_TIME}/wildtype_000000/kb"
VARIANT_METADATA_DIR="out/${SUBMISSION_TIME}/wildtype_000000/metadata"


##### Create output directory #####

mkdir -p "${SIM_OUT_DATA_DIR}"


##### Create metadata directory #####

mkdir -p "${METADATA_DIR}"


##### Create directories to match multi-job file structure #####

mkdir -p "${VARIANT_KB_DIR}"
mkdir -p "${VARIANT_METADATA_DIR}"


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


##### Create knowledgebases #####

make justKb FIXTURES_KBDIR="${KB_DIR}"


##### Run simulation #####

PYTHONPATH="$PWD:$PYTHONPATH" WC_KBLOCATION="\"${KB_FIT}\"" python2.7 runscripts/runSimulation.py "${SUBMISSION_TIME}"

# If the simulation didn't complete successfully, don't run analysis
if [ "$?" -ne "0" ]; then
	exit 1
fi

##### Single simulation analysis #####

SINGLE_ANALYSIS_SCRIPTS=$(PYTHONPATH="$PWD:$PYTHONPATH" python2.7 -c "from models.ecoli.sim.simulation import EcoliSimulation; EcoliSimulation.printAnalysisSingleFiles()")

PLOT_OUT_DATA_DIR="out/${SUBMISSION_TIME}/wildtype_000000/${SEED_DIR}/plotOut"

mkdir -p $PLOT_OUT_DATA_DIR

echo "+++++ Processing $SIM_OUT_DATA_DIR +++++"

SAVEIFS=$IFS
IFS=$(echo -en "\n\b")
for SINGLE_ANALYSIS_SCRIPT in $SINGLE_ANALYSIS_SCRIPTS; do
	if [ "$(basename $SINGLE_ANALYSIS_SCRIPT)" = "__init__.py" ]; then
		continue
	fi

	OUT_NAME=$(basename $SINGLE_ANALYSIS_SCRIPT | sed 's/\.py//g')

	echo "Running $(basename $SINGLE_ANALYSIS_SCRIPT)"

	PYTHONPATH="$PWD:$PYTHONPATH" python2.7 "$SINGLE_ANALYSIS_SCRIPT" "$SIM_OUT_DATA_DIR" "$PLOT_OUT_DATA_DIR" "${OUT_NAME}" --kbFile "${KB_FIT}"
done
IFS=$SAVEIFS

##### Compress kb fixtures to save space #####

echo "+++++ Compressing KB files +++++"

KB_FILES=$(find ${KB_DIR} -type f | sort)

# TODO: Decide if you want to do this (might be good to leave the links for documentation?)
# Remove symlinks (their targets won't exist anymore)
#find ${KB_DIR} -type l -exec unlink {} \;

for KB_FILE in $KB_FILES; do
	bzip2 "${KB_FILE}"
done

echo
