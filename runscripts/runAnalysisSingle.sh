#!/bin/bash

if [ "$#" -ne 3 ]; then
	echo "Usage: $0 SIM_OUT_DIR PLOT_OUT_DIR SCRIPTS_DIR" >& 2
	exit 1
fi

SIM_OUT_DIR="$1"
PLOT_OUT_DIR="${2%/}/"
SCRIPTS_DIR="$3"

SIM_OUT_DATA_DIRS=$(find $SIM_OUT_DIR -name "*.hdf" -print0 | xargs -0 -I {} dirname {} | uniq | sort)

SCRIPTS=$(find $SCRIPTS_DIR -name "*.py" | sort)

for SIM_OUT_DATA_DIR in $SIM_OUT_DATA_DIRS; do
	PLOT_OUT_DATA_DIR=${SIM_OUT_DATA_DIR/$SIM_OUT_DIR/$PLOT_OUT_DIR}

	mkdir -p $PLOT_OUT_DATA_DIR
	
	echo "+++++ Processing $SIM_OUT_DATA_DIR +++++"
	
	for SCRIPT in $SCRIPTS; do
		if [ "$(basename $SCRIPT)" = "__init__.py" ]; then
			continue
		fi

		OUT_NAME=$(basename $SCRIPT | sed 's/.py//g')

		echo "Running $(basename $SCRIPT)"

		PYTHONPATH="$PWD:$PYTHONPATH" python2.7 $SCRIPT $SIM_OUT_DATA_DIR $PLOT_OUT_DATA_DIR ${OUT_NAME}.pdf
	done

	echo
done