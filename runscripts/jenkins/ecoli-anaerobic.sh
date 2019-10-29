HOST=$1
NAME=$2
PORT=$3
PASSWORD=$4

set -e

module load wcEcoli/sherlock2
pyenv local wcEcoli-paper

make clean
make compile

sh runscripts/jenkins/fireworks-config.sh $HOST $NAME $PORT $PASSWORD

echo y | lpad reset

PYTHONPATH=$PWD DESC="Anaerobic." VARIANT="condition" FIRST_VARIANT_INDEX=1 LAST_VARIANT_INDEX=1 SINGLE_DAUGHTERS=1 N_GENS=8 MASS_DISTRIBUTION=0 COMPRESS_OUTPUT=1 python runscripts/fw_queue.py

while [ $(lpad get_fws -s READY -d count) -ge 1 ]; do
  PYTHONPATH=$PWD rlaunch singleshot
done

N_FAILS=$(lpad get_fws -s FIZZLED -d count)

test $N_FAILS = 0

mv out/2* /scratch/PI/mcovert/wc_ecoli/anaerobic/
