set -e

### ---------------------------------------------------------------------------
### Edit the pyenv in setup-environment.sh to make the PR build temporarily use
### another pyenv for testing (eg. wcEcoli3-staging). Revert the change before
### merging the PR into master to prevent changing it for other Jenkins builds.
### ---------------------------------------------------------------------------
source runscripts/jenkins/setup-environment.sh
sh runscripts/jenkins/fireworks-config.sh "pr$EXECUTOR_NUMBER"

echo y | lpad reset

DESC="2 generations completion test." OPERONS=off WC_ANALYZE_FAST=1 SINGLE_DAUGHTERS=1 N_GENS=2 MASS_DISTRIBUTION=0 \
	PARALLEL_PARCA=1 COMPRESS_OUTPUT=0 PLOTS=ACTIVE BUILD_CAUSALITY_NETWORK=1 RAISE_ON_TIME_LIMIT=1 python runscripts/fireworks/fw_queue.py

bash runscripts/jenkins/run-fireworks.sh

runscripts/jenkins/runscript-checks.sh

rm -fr out/*
