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

CODE_DIR="$PBS_O_WORKDIR" # Assumes job submission from wcEcoli

if [ -z "${VARIANT}" ]; then
	KB_DIR="${CODE_DIR}/out/${SUBMISSION_TIME}/kb"
else
	if [ -z "${VARIANT_ID}" ]; then
		echo "VARIANT_ID environmental variable must be set" >&2
		exit 1
	fi
	KB_DIR="${CODE_DIR}/out/${SUBMISSION_TIME}/${VARIANT}_$(printf "%06d" ${VARIANT_ID})/kb"
fi


echo KB_DIR ${KB_DIR}

runprogram()
{
	echo "Compressing files"

	KB_FILES=$(find ${KB_DIR} -type f | sort)

	# TODO: Decide if you want to do this (might be good to leave the links for documentation?)
	# Remove symlinks (their targets won't exist anymore)
	#find ${KB_DIR} -type l -exec unlink {} \;

	for KB_FILE in $KB_FILES; do
		bzip2 "${KB_FILE}"
	done

}

early()
{
	echo
	echo "##### WARNING: EARLY TERMINATION #####"
	echo
}

trap "early; stageout" 2 9 15

runprogram
