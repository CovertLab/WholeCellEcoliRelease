HOST=$1
NAME=$2
PORT=$3
PASSWORD=$4

set -e

module load wcEcoli/python3

### -------------------------------------------------------------------
### Edit this line to make the PR build use another pyenv like wcEcoli3-staging.
### Revert it to `wcEcoli3` before merging the PR into master.
### -------------------------------------------------------------------
pyenv local wcEcoli3

make clean compile

# Get mypy type checker warnings now but defer failing on its error detections.
set +e
runscripts/debug/mypy.sh
MYPY_FAILED=$?
set -e

PYTHONPATH=$PWD:$PYTHONPATH pytest --cov=wholecell --cov-report xml \
    --junitxml=unittests.xml

sh runscripts/jenkins/fireworks-config.sh $HOST $NAME $PORT $PASSWORD

echo y | lpad reset

PYTHONPATH=$PWD DESC="2 generations completion test." WC_ANALYZE_FAST=1 SINGLE_DAUGHTERS=1 N_GENS=2 MASS_DISTRIBUTION=0 \
	PARALLEL_PARCA=1 COMPRESS_OUTPUT=0 PLOTS=ACTIVE BUILD_CAUSALITY_NETWORK=1 RAISE_ON_TIME_LIMIT=1 python runscripts/fireworks/fw_queue.py

PYTHONPATH=$PWD rlaunch rapidfire --nlaunches 0

N_FAILS=$(lpad get_fws -s FIZZLED -d count)

if [ $N_FAILS -gt 0 ]; then
  mv out/2* /scratch/PI/mcovert/wc_ecoli/failed/
fi

git status | head -1

test $N_FAILS = 0 -a $MYPY_FAILED = 0

rm -fr out/*
