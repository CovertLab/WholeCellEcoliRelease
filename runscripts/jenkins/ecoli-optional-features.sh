set -e

source runscripts/jenkins/setup-environment.sh
sh runscripts/jenkins/fireworks-config.sh optional

echo y | lpad reset

DESC="No tRNA Charging" TRNA_CHARGING=0 N_GENS=8 \
  PARALLEL_PARCA=1 SINGLE_DAUGHTERS=1 COMPRESS_OUTPUT=1 RAISE_ON_TIME_LIMIT=1 \
  PLOTS=ACTIVE WC_ANALYZE_FAST=1 \
  python runscripts/fireworks/fw_queue.py
# Better regulation of replisome subunits should allow MECHANISTIC_REPLISOME=1
DESC="ppGpp regulation" PPGPP_REGULATION=1 MECHANISTIC_REPLISOME=0 N_GENS=8 \
  PARALLEL_PARCA=1 SINGLE_DAUGHTERS=1 COMPRESS_OUTPUT=1 RAISE_ON_TIME_LIMIT=1 \
  PLOTS=ACTIVE WC_ANALYZE_FAST=1 \
  python runscripts/fireworks/fw_queue.py
DESC="Superhelical Densities" SUPERHELICAL_DENSITIES=1 N_GENS=8 \
  PARALLEL_PARCA=1 SINGLE_DAUGHTERS=1 COMPRESS_OUTPUT=1 RAISE_ON_TIME_LIMIT=1 \
  PLOTS=ACTIVE WC_ANALYZE_FAST=1 \
  python runscripts/fireworks/fw_queue.py
DESC="Causality Network" BUILD_CAUSALITY_NETWORK=1 N_GENS=2 SEED=$RANDOM \
  PARALLEL_PARCA=1 SINGLE_DAUGHTERS=1 COMPRESS_OUTPUT=1 RAISE_ON_TIME_LIMIT=1 \
  WC_ANALYZE_FAST=1 \
  python runscripts/fireworks/fw_queue.py

bash runscripts/jenkins/run-fireworks.sh

runscripts/jenkins/save_output.sh out/ /scratch/PI/mcovert/wc_ecoli/optional_features/
