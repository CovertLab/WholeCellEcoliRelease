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
KB_FIT="${KB_DIR}/KnowledgeBase_Fit.cPickle"

RESULTS_DIR="${CODE_DIR}/out/${SUBMISSION_TIME}"

# mkdir -p "$RESULTS_DIR"

SEED=$(printf "%06d" $(($ARRAY_ID - 1)))
OUTPUT_LOG_BASE_NAME="simShellLog"
OUTPUT_LOG_FILE="${PBS_O_WORKDIR}/${OUTPUT_LOG_BASE_NAME}.${SUBMISSION_TIME}:$SEED"

echo WORK_DIR $WORK_DIR
echo CODE_DIR $CODE_DIR
echo OUTPUT_LOG_FILE $OUTPUT_LOG_FILE

stagein()
{
	echo
	echo "Copying files to work directory ${WORK_DIR}"

	cd ${WORK_DIR}
	scp -r ${KBECOLI_DIR} .

	mkdir $(basename $CODE_DIR)
	cd $(basename $CODE_DIR)
	scp -r ${KB_FIT} .
	# scp -r ${CODE_DIR}/fixtures .
	scp -r ${CODE_DIR}/runscripts .
	scp -r ${CODE_DIR}/user .
	scp -r ${CODE_DIR}/wholecell .

	mkdir -p "out/${SUBMISSION_TIME}/${SEED}/simOut"

}

runprogram()
{
	echo "Running"

	cd ${WORK_DIR}/$(basename $CODE_DIR)
	WC_SEED=${SEED} \
	WC_STATES=${WC_STATES} \
	WC_PROCESSES=${WC_PROCESSES} \
	WC_LISTENERS=${WC_LISTENERS} \
	WC_HOOKS=${WC_HOOKS} \
	WC_LENGTHSEC=${WC_LENGTHSEC} \
	WC_TIMESTEPSEC=${WC_TIMESTEPSEC} \
	WC_LOGTOSHELL=${WC_LOGTOSHELL} \
	WC_SHELLCOLUMNSHEADERS=${WC_SHELLCOLUMNSHEADERS} \
	WC_LOGTODISKEVERY=${WC_LOGTODISKEVERY} \
	WC_KBLOCATION="\"${WORK_DIR}/$(basename $CODE_DIR)/KnowledgeBase_Fit.cPickle\"" python2.7 runscripts/runSimulationJob.py "${SUBMISSION_TIME}" 2>&1 | tee -a "${OUTPUT_LOG_FILE}"
}

stageout()
{
	echo "Transferring files back"

	cd ${WORK_DIR}/$(basename $CODE_DIR)
	mkdir -p "$RESULTS_DIR/${SEED}/simOut"
	scp -r "out/${SUBMISSION_TIME}/${SEED}/simOut" "$RESULTS_DIR/${SEED}"
	mv "${OUTPUT_LOG_FILE}" "${RESULTS_DIR}/${SEED}/simOut/${OUTPUT_LOG_BASE_NAME}"

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

exit 0