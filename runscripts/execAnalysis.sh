#!/bin/bash

if [ "$#" -ne 3 ]; then
	echo "Usage: $0 MODEL_LEVEL KB_DIR SIM_DIR" >& 2
	exit 1
fi

MODEL_LEVEL="$1"
KB_DIR="$2"
SIM_DIR="$3"

SIM_OUT_DATA_DIR="${SIM_DIR}/model_level_${MODEL_LEVEL}/simOut"
METADATA_DIR="${SIM_DIR}/model_level_${MODEL_LEVEL}/metadata"
KB_FILE="${KB_DIR}/KnowledgeBase_Most_Fit.cPickle"

SINGLE_ANALYSIS_SCRIPTS="$(cat ${METADATA_DIR}/singleAnalysis.list)"
PLOT_OUT_DATA_DIR="${SIM_DIR}/model_level_${MODEL_LEVEL}/plotOut"

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
