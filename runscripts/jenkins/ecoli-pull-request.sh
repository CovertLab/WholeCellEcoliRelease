HOST=$1
NAME=$2
PORT=$3
PASSWORD=$4

set -e

module load wcEcoli/sherlock2

### -------------------------------------------------------------------
### Edit this line to make the PR build use a new pyenv.
### Revert to `pyenv local wcEcoli2` before merging the PR into master.
### -------------------------------------------------------------------
pyenv local wcEcoli2

make clean
make compile

PYTHONPATH=$PWD:$PYTHONPATH nosetests --with-xunit --with-coverage --cover-package=wholecell --cover-xml

sh runscripts/jenkins/fireworks-config.sh $HOST $NAME $PORT $PASSWORD

echo y | lpad reset

PYTHONPATH=$PWD DESC="2 generations completion test." WC_ANALYZE_FAST=1 SINGLE_DAUGHTERS=1 N_GENS=2 MASS_DISTRIBUTION=0 \
	PARALLEL_PARCA=1 COMPRESS_OUTPUT=0 PLOTS=ACTIVE BUILD_CAUSALITY_NETWORK=1 RAISE_ON_TIME_LIMIT=1 python runscripts/fireworks/fw_queue.py

PYTHONPATH=$PWD rlaunch rapidfire --nlaunches 0

N_FAILS=$(lpad get_fws -s FIZZLED -d count)

if [ $N_FAILS -gt 0 ]; then
  mv out/2* /scratch/PI/mcovert/wc_ecoli/failed/
fi

test $N_FAILS = 0

rm -fr out/*
