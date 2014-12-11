# Merge stdin and stderr
#PBS -j oe

# Plenty of time to run
#PBS -l walltime=24:00:00

# Plenty of memory
#PBS -l mem=2G


echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo PBS: PYTHONPATH = $PYTHONPATH
echo ------------------------------------------------------

if [ -z "$SUBMISSION_TIME" ]; then
	echo "SUBMISSION_TIME environmental variable must be set" >&2
	exit 1
fi

if [ -z "$ARRAY_ID" ]; then
	echo "ARRAY_ID environmental variable must be set" >&2
	exit 1
fi

if [ -d "/state/partition1" ]; then
	WORK_DIR="/state/partition1"
else
	WORK_DIR="/tmp"
fi

WORK_DIR="${WORK_DIR}/${SUBMISSION_TIME}.${PBS_JOBID}.${ARRAY_ID}"

mkdir -p "$WORK_DIR"

CODE_DIR="$PBS_O_WORKDIR" # Assumes job submission from wcEcoli
KBECOLI_DIR="${CODE_DIR}/../kbEcoli"
KB_DIR="${CODE_DIR}/out/${SUBMISSION_TIME}/kb"
KB_FIT="${KB_DIR}/KnowledgeBase_Most_Fit.cPickle"

SEED=$(printf "%06d" $(($ARRAY_ID - 1)))

RESULTS_DIR="${CODE_DIR}/out/${SUBMISSION_TIME}"
SPECIFIC_RESULTS_DIR="${RESULTS_DIR}/${SEED}/simOut"

# PLOTS_DIR="${CODE_DIR}/out/plotOut"
# SUBMISSION_PLOTS_DIR="${PLOTS_DIR}/${SUBMISSION_TIME}"
SPECIFIC_PLOTS_DIR="${RESULTS_DIR}/${SEED}/plotOut"

mkdir -p "$SPECIFIC_PLOTS_DIR"

MY_SPECIFIC_RESULTS_DIR="${WORK_DIR}/$(basename $CODE_DIR)/out/${SUBMISSION_TIME}/${SEED}/simOut"
MY_SPECIFIC_PLOTS_DIR="${WORK_DIR}/$(basename $CODE_DIR)/out/${SUBMISSION_TIME}/${SEED}/plotOut"


# MY_SPECIFIC_RESULTS_DIR="${WORK_DIR}/$(basename $CODE_DIR)/out/simOut/${SUBMISSION_TIME}/$(basename $SPECIFIC_RESULTS_DIR)"
# MY_SPECIFIC_PLOTS_DIR="${WORK_DIR}/$(basename $CODE_DIR)/out/plotOut/${SUBMISSION_TIME}/$(basename $SPECIFIC_PLOTS_DIR)"


OUTPUT_LOG_BASE_NAME="analysisSingleLog"
OUTPUT_LOG_FILE="${PBS_O_WORKDIR}/${OUTPUT_LOG_BASE_NAME}.${SUBMISSION_TIME}:$SEED"

echo WORK_DIR $WORK_DIR
echo CODE_DIR $CODE_DIR
echo SPECIFIC_RESULTS_DIR $SPECIFIC_RESULTS_DIR
echo SPECIFIC_PLOTS_DIR $SPECIFIC_PLOTS_DIR
echo MY_SPECIFIC_RESULTS_DIR $MY_SPECIFIC_RESULTS_DIR
echo MY_SPECIFIC_PLOTS_DIR $MY_SPECIFIC_PLOTS_DIR
echo OUTPUT_LOG_FILE $OUTPUT_LOG_FILE

stagein()
{
	echo
	echo "Copying files to work directory ${WORK_DIR}"

	cd ${WORK_DIR}

	scp -r ${KBECOLI_DIR} .

	mkdir $(basename $CODE_DIR)
	cd $(basename $CODE_DIR)
	scp -r ${CODE_DIR}/runscripts .
	scp -r ${CODE_DIR}/user .
	scp -r ${CODE_DIR}/wholecell .
	scp -r ${CODE_DIR}/models .
	scp -r ${CODE_DIR}/reconstruction .

	mkdir -p "out/${SUBMISSION_TIME}/${SEED}"

	cd "out/${SUBMISSION_TIME}"
	scp -r "${KB_DIR}" .

	cd "${SEED}"
	scp -r ${SPECIFIC_RESULTS_DIR} .
}

runprogram()
{
	echo "Running"

	mkdir -p "$MY_SPECIFIC_PLOTS_DIR"

	cd ${WORK_DIR}/$(basename $CODE_DIR)
	SCRIPTS=$(PYTHONPATH="$PWD:$PYTHONPATH" python2.7 -c "from models.ecoli.sim.simulation import EcoliSimulation; EcoliSimulation.printAnalysisSingleFiles()")

	SAVEIFS=$IFS
	IFS=$(echo -en "\n\b")
	for SCRIPT in $SCRIPTS; do
		if [ "$(basename $SCRIPT)" = "__init__.py" ]; then
			continue
		fi

		OUT_NAME=$(basename $SCRIPT | sed 's/.py//g')

		echo "Running $(basename $SCRIPT)"

		PYTHONPATH="${WORK_DIR}/$(basename $CODE_DIR):$PYTHONPATH" python2.7 "$SCRIPT" "$MY_SPECIFIC_RESULTS_DIR" "$MY_SPECIFIC_PLOTS_DIR" "${OUT_NAME}.svg" --kbFile "${WORK_DIR}/$(basename $CODE_DIR)/out/${SUBMISSION_TIME}/kb/KnowledgeBase_Most_Fit.cPickle"
	done 2>&1 | tee -a "${OUTPUT_LOG_FILE}"
	IFS=$SAVEIFS
}

stageout()
{
	echo "Transferring files back"

	scp -r $MY_SPECIFIC_PLOTS_DIR "${RESULTS_DIR}/${SEED}"
	mv "${OUTPUT_LOG_FILE}" "${RESULTS_DIR}/${SEED}/plotOut/${OUTPUT_LOG_BASE_NAME}"

	echo "Cleaning up"
	cd /
	rm -fr "${WORK_DIR}"
}

early()
{
	echo
	echo "##### WARNING: EARLY TERMINATION #####"
	echo
}

trap "early; stageout" 2 9 15

stagein
runprogram
stageout
