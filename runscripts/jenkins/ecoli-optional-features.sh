HOST=$1
NAME=$2
PORT=$3
PASSWORD=$4

set -e

module load wcEcoli/sherlock2
pyenv local wcEcoli2

make clean
make compile

sh runscripts/jenkins/fireworks-config.sh $HOST $NAME $PORT $PASSWORD

echo y | lpad reset

PYTHONPATH=$PWD DESC="No tRNA Charging" SINGLE_DAUGHTERS=1 N_GENS=8 TRNA_CHARGING=0 COMPRESS_OUTPUT=1 PLOTS=ACTIVE RAISE_ON_TIME_LIMIT=1 python runscripts/fireworks/fw_queue.py
PYTHONPATH=$PWD DESC="ppGpp regulation" PPGPP_REGULATION=1 SINGLE_DAUGHTERS=1 N_GENS=8 COMPRESS_OUTPUT=1 PLOTS=ACTIVE RAISE_ON_TIME_LIMIT=1 python runscripts/fireworks/fw_queue.py
PYTHONPATH=$PWD DESC="Causality Network" SINGLE_DAUGHTERS=1 N_GENS=2 BUILD_CAUSALITY_NETWORK=1 COMPRESS_OUTPUT=1 RAISE_ON_TIME_LIMIT=1 python runscripts/fireworks/fw_queue.py

PYTHONPATH=$PWD rlaunch rapidfire --nlaunches 0

N_FAILS=$(lpad get_fws -s FIZZLED -d count)

if [ $N_FAILS -gt 0 ]; then
  mv out/2* /scratch/PI/mcovert/wc_ecoli/failed/
fi

test $N_FAILS = 0

mv out/2* /scratch/PI/mcovert/wc_ecoli/optional_features/
