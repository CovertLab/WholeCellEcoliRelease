#!/bin/bash

trap "exit 130" 2 9 15

SUBMISSION_TIME=$(date "+%Y%m%d.%H%M%S.%N")
WC_SEED=${WC_SEED:-0}
SEED_DIR=$(printf "%06d" $(($WC_SEED)))
SIM_OUT_DATA_DIR="out/simOut/${SUBMISSION_TIME}/${SEED_DIR}"


##### Create output directory #####
mkdir -p "${SIM_OUT_DATA_DIR}"


##### Create knowledgebases (unfit and fit) #####
python2.7 runscripts/createKbs.py --outputDirectory "out/simOut/${SUBMISSION_TIME}/kb"


##### Run simulation #####
python2.7 runscripts/runSimulation.py "${SUBMISSION_TIME}"

# If the simulation didn't complete successfully, don't run analysis
if [ "$?" -ne "0" ]; then
	exit 1
fi

##### Single simulation analysis #####

SINGLE_ANALYSIS_SCRIPTS_DIR="wholecell/analysis/single"
SINGLE_ANALYSIS_SCRIPTS=$(find $SINGLE_ANALYSIS_SCRIPTS_DIR -name "*.py" | sort)

PLOT_OUT_DATA_DIR="out/plotOut/${SUBMISSION_TIME}/${SEED_DIR}"

mkdir -p $PLOT_OUT_DATA_DIR

echo "+++++ Processing $SIM_OUT_DATA_DIR +++++"

for SINGLE_ANALYSIS_SCRIPT in $SINGLE_ANALYSIS_SCRIPTS; do
	if [ "$(basename $SINGLE_ANALYSIS_SCRIPT)" = "__init__.py" ]; then
		continue
	fi

	OUT_NAME=$(basename $SINGLE_ANALYSIS_SCRIPT | sed 's/.py//g')

	echo "Running $(basename $SINGLE_ANALYSIS_SCRIPT)"

	python2.7 $SINGLE_ANALYSIS_SCRIPT $SIM_OUT_DATA_DIR $PLOT_OUT_DATA_DIR ${OUT_NAME}.pdf
done

echo
