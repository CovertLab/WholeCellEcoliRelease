HOST=$1
NAME=$2
PORT=$3
PASSWORD=$4

set -e

runscripts/jenkins/purge.sh anaerobic 10

module load wcEcoli/python3
pyenv local wcEcoli3

make clean
make compile

sh runscripts/jenkins/fireworks-config.sh $HOST $NAME $PORT $PASSWORD

echo y | lpad reset

PYTHONPATH=$PWD DESC="Anaerobic." VARIANT="condition" FIRST_VARIANT_INDEX=1 LAST_VARIANT_INDEX=1 SINGLE_DAUGHTERS=1 N_GENS=8 MASS_DISTRIBUTION=0 COMPRESS_OUTPUT=1 PLOTS=ACTIVE RAISE_ON_TIME_LIMIT=1 python runscripts/fireworks/fw_queue.py

while [ $(lpad get_fws -s READY -d count) -ge 1 ]; do
  PYTHONPATH=$PWD rlaunch singleshot
done

N_FAILS=$(lpad get_fws -s FIZZLED -d count)

if [ $N_FAILS -gt 0 ]; then
  mv out/2* /scratch/PI/mcovert/wc_ecoli/failed/
fi

test $N_FAILS = 0

runscripts/jenkins/save_output.sh out/ /scratch/PI/mcovert/wc_ecoli/anaerobic/
