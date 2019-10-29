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

PYTHONPATH=$PWD DESC="With AAs." VARIANT="condition" FIRST_VARIANT_INDEX=2 LAST_VARIANT_INDEX=2 SINGLE_DAUGHTERS=1 N_GENS=8 MASS_DISTRIBUTION=0 COMPRESS_OUTPUT=1 python runscripts/fw_queue.py

PYTHONPATH=$PWD rlaunch rapidfire --nlaunches 0

N_FAILS=$(lpad get_fws -s FIZZLED -d count)

test $N_FAILS = 0

mv out/2* /scratch/PI/mcovert/wc_ecoli/with_aa/
