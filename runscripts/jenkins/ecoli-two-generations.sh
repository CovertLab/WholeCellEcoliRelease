HOST=$1
NAME=$2
PORT=$3
PASSWORD=$4

set -e

source runscripts/jenkins/setup-environment.sh
sh runscripts/jenkins/fireworks-config.sh $HOST $NAME $PORT $PASSWORD

echo y | lpad reset

DESC="2 generations completion test." SINGLE_DAUGHTERS=1 N_GENS=2 MASS_DISTRIBUTION=0 COMPRESS_OUTPUT=1 RAISE_ON_TIME_LIMIT=1 python runscripts/fireworks/fw_queue.py

bash runscripts/jenkins/run-fireworks.sh

rm -fr out/*

# Run everything with WC_ANALYZE_FAST and PARALLEL_PARCA

echo y | lpad reset

DESC="2 generations completion test." WC_ANALYZE_FAST=1 SINGLE_DAUGHTERS=1 N_GENS=2 PARALLEL_PARCA=1 PLOTS=ACTIVE RAISE_ON_TIME_LIMIT=1 python runscripts/fireworks/fw_queue.py

bash runscripts/jenkins/run-fireworks.sh

rm -fr out/*
