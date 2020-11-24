HOST=$1
NAME=$2
PORT=$3
PASSWORD=$4

set -e

### ---------------------------------------------------------------------------
### Edit the pyenv in setup-environment.sh to make the PR build temporarily use
### another pyenv for testing (eg. wcEcoli3-staging). Revert the change before
### merging the PR into master to prevent changing it for other Jenkins builds.
### ---------------------------------------------------------------------------
source runscripts/jenkins/setup-environment.sh
sh runscripts/jenkins/fireworks-config.sh $HOST $NAME $PORT $PASSWORD

echo y | lpad reset

DESC="2 generations completion test." WC_ANALYZE_FAST=1 SINGLE_DAUGHTERS=1 N_GENS=2 MASS_DISTRIBUTION=0 \
	PARALLEL_PARCA=1 COMPRESS_OUTPUT=0 PLOTS=ACTIVE BUILD_CAUSALITY_NETWORK=1 RAISE_ON_TIME_LIMIT=1 python runscripts/fireworks/fw_queue.py

rlaunch rapidfire --nlaunches 0

N_FAILS=$(lpad get_fws -s FIZZLED -d count)

if [ $N_FAILS -gt 0 ]; then
  lpad get_fws -s FIZZLED
  mv out/2* /scratch/PI/mcovert/wc_ecoli/failed/
fi

git status | head -1

test $N_FAILS = 0

rm -fr out/*
