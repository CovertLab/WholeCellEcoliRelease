set -e

source runscripts/jenkins/setup-environment.sh
sh runscripts/jenkins/fireworks-config.sh anaerobic

echo y | lpad reset

DESC="Anaerobic." VARIANT="condition" FIRST_VARIANT_INDEX=1 LAST_VARIANT_INDEX=1 SINGLE_DAUGHTERS=1 N_GENS=8 MASS_DISTRIBUTION=0 COMPRESS_OUTPUT=1 PLOTS=ACTIVE RAISE_ON_TIME_LIMIT=1 python runscripts/fireworks/fw_queue.py

bash runscripts/jenkins/run-fireworks.sh

runscripts/jenkins/save_output.sh out/ /scratch/PI/mcovert/wc_ecoli/anaerobic/
