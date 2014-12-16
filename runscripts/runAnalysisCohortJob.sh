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

if [ -d "/state/partition1" ]; then
	WORK_DIR="/state/partition1"
else
	WORK_DIR="/tmp"
fi

WORK_DIR="${WORK_DIR}/${SUBMISSION_TIME}.${PBS_JOBID}"

mkdir -p "$WORK_DIR"

CODE_DIR="$PBS_O_WORKDIR" # Assumes job submission from wcEcoli
RESULTS_DIR="${CODE_DIR}/out/simOut"
SUBMISSION_RESULTS_DIR="${RESULTS_DIR}/${SUBMISSION_TIME}"

PLOTS_DIR="${CODE_DIR}/out/plotOut"
SUBMISSION_PLOTS_DIR="${PLOTS_DIR}/${SUBMISSION_TIME}"
COHORT_PLOTS_DIR="${SUBMISSION_PLOTS_DIR}/cohort"

MY_SUBMISSION_RESULTS_DIR="${WORK_DIR}/$(basename $CODE_DIR)/out/simOut/${SUBMISSION_TIME}"
MY_SUBMISSION_PLOTS_DIR="${WORK_DIR}/$(basename $CODE_DIR)/out/plotOut/${SUBMISSION_TIME}"

OUTPUT_LOG_BASE_NAME="analysisCohortLog"
OUTPUT_LOG_FILE="${PBS_O_WORKDIR}/${OUTPUT_LOG_BASE_NAME}.${SUBMISSION_TIME}"

echo WORK_DIR $WORK_DIR
echo CODE_DIR $CODE_DIR
echo MY_SUBMISSION_RESULTS_DIR $MY_SUBMISSION_RESULTS_DIR
echo MY_SUBMISSION_PLOTS_DIR $MY_SUBMISSION_PLOTS_DIR
echo OUTPUT_LOG_FILE $OUTPUT_LOG_FILE

stagein()
{
	echo
	echo "Copying files to work directory ${WORK_DIR}"

	cd ${WORK_DIR}

	mkdir $(basename $CODE_DIR)
	cd $(basename $CODE_DIR)
	scp -r ${CODE_DIR}/fixtures .
	scp -r ${CODE_DIR}/runscripts .
	scp -r ${CODE_DIR}/user .
	scp -r ${CODE_DIR}/wholecell .
	scp -r ${CODE_DIR}/models .
	scp -r ${CODE_DIR}/reconstruction .

	mkdir -p out/simOut/${SUBMISSION_TIME}
	scp -r ${SUBMISSION_RESULTS_DIR} out/simOut
}

runprogram()
{
	echo "Running"

	mkdir -p "$MY_SUBMISSION_PLOTS_DIR"

	cd ${WORK_DIR}/$(basename $CODE_DIR)
	SCRIPTS=$(PYTHONPATH="$PWD:$PYTHONPATH" python2.7 -c "from models.ecoli.sim.simulation import EcoliSimulation; EcoliSimulation.printAnalysisCohortFiles()")

	for SCRIPT in $SCRIPTS; do
		if [ "$(basename $SCRIPT)" = "__init__.py" ]; then
			continue
		fi

		OUT_NAME=$(basename $SCRIPT | sed 's/.py//g')

		echo "Running $(basename $SCRIPT)"

		python2.7 $SCRIPT $MY_SUBMISSION_RESULTS_DIR $MY_SUBMISSION_PLOTS_DIR ${OUT_NAME}.svg
	done 2>&1 | tee -a "${OUTPUT_LOG_FILE}"
}

stageout()
{
	echo "Transferring files back"

	scp -r "${MY_SUBMISSION_PLOTS_DIR}" "${COHORT_PLOTS_DIR}"
	mv "${OUTPUT_LOG_FILE}" "${SUBMISSION_PLOTS_DIR}/${OUTPUT_LOG_BASE_NAME}"

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
