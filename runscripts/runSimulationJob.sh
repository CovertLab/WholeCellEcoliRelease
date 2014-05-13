# Merge stdin and stderr
#PBS -j oe

# Plenty of time to run
#PBS -l walltime=24:00:00

# Plenty of memory
#PBS -l mem=2G

# Job array IDs
# #PBS -t 0-1


sleep 10
date

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
echo PBS: job array id = $PBS_ARRAYID
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

WORK_DIR="${WORK_DIR}/${SUBMISSION_TIME}"

if [ -z "$PBS_ARRAYID" ]; then
	echo "PBS_ARRAYID not set"
else
	WORK_DIR="${WORK_DIR}/${PBS_ARRAYID}"
	echo "PBS_ARRAYID is set. WORK_DIR is $WORK_DIR."
fi

mkdir -p "$WORK_DIR"

CODE_DIR="$PBS_O_WORKDIR" # Assumes job submission from wcEcoli
KB_DIR="${CODE_DIR}/../kbEcoli"

RESULTS_DIR="${CODE_DIR}/out/simOut"

mkdir -p "$RESULTS_DIR"

echo WORK_DIR $WORK_DIR
echo CODE_DIR $CODE_DIR

stagein()
{
	echo
	echo "Copying files to work directory ${WORK_DIR}"

	cd ${WORK_DIR}
	scp -r ${CODE_DIR} .
	scp -r ${KB_DIR} .
}

runprogram()
{
	echo "Running"

	cd ${WORK_DIR}/$(basename $CODE_DIR)
	python2.7 runscripts/runSimulationJob.py "${SUBMISSION_TIME}"
}

stageout()
{
	echo "Transferring files back"

	cd ${WORK_DIR}/$(basename $CODE_DIR)
	scp -r "out/simOut/${SUBMISSION_TIME}" "$RESULTS_DIR"

	cd /
	rm -fr "${WORK_DIR}/$(basename $CODE_DIR)"
	rm -fr "${WORK_DIR}/kbEcoli"
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
