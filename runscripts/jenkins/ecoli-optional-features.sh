HOST=$1
NAME=$2
PORT=$3
PASSWORD=$4

set -e

module load wcEcoli/sherlock2
pyenv local wcEcoli2
export PYTHONPATH=$PWD

make clean
make compile

sh runscripts/jenkins/fireworks-config.sh $HOST $NAME $PORT $PASSWORD

echo y | lpad reset

DESC="No tRNA Charging" TRNA_CHARGING=0 N_GENS=8 \
  PARALLEL_PARCA=1 SINGLE_DAUGHTERS=1 COMPRESS_OUTPUT=1 RAISE_ON_TIME_LIMIT=1 \
  PLOTS=ACTIVE WC_ANALYZE_FAST=1 \
  python runscripts/fireworks/fw_queue.py
DESC="ppGpp regulation" PPGPP_REGULATION=1 N_GENS=8 \
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

# Commented rapidfire command below produces seg fault after 2 hr and 10 min (see #764)
# Could replace singleshot loop with rapidfire if fixed
# Singleshot might seg fault as well for long single tasks over 2 hr and 10 min

# rlaunch rapidfire --nlaunches 0
while [ $(lpad get_fws -s READY -d count) -ge 1 ]; do
  rlaunch singleshot
done

N_FAILS=$(lpad get_fws -s FIZZLED -d count)

if [ $N_FAILS -gt 0 ]; then
  mv out/2* /scratch/PI/mcovert/wc_ecoli/failed/
fi

test $N_FAILS = 0

runscripts/jenkins/save_output.sh out/ /scratch/PI/mcovert/wc_ecoli/optional_features/
