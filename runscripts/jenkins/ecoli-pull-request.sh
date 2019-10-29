HOST=$1
NAME=$2
PORT=$3
PASSWORD=$4

set -e

module load wcEcoli/sherlock2
pyenv local wcEcoli-paper

make clean
make compile

PYTHONPATH=$PWD:$PYTHONPATH nosetests --with-xunit --with-coverage --cover-package=wholecell --cover-xml

sh runscripts/jenkins/fireworks-config.sh $HOST $NAME $PORT $PASSWORD

echo y | lpad reset

PYTHONPATH=$PWD DESC="2 generations completion test." WC_ANALYZE_FAST=1 SINGLE_DAUGHTERS=1 N_GENS=2 MASS_DISTRIBUTION=0 PARALLEL_FITTER=1 COMPRESS_OUTPUT=0 python runscripts/fw_queue.py

PYTHONPATH=$PWD rlaunch rapidfire --nlaunches 0

N_FAILS=$(lpad get_fws -s FIZZLED -d count)

test $N_FAILS = 0
