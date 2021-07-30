set -e

### ---------------------------------------------------------------------------
### Edit the pyenv in setup-environment.sh to make the PR build temporarily use
### another pyenv for testing (eg. wcEcoli3-staging). Revert the change before
### merging the PR into master to prevent changing it for other Jenkins builds.
### ---------------------------------------------------------------------------
source runscripts/jenkins/setup-environment.sh
sh runscripts/jenkins/fireworks-config.sh "pr$EXECUTOR_NUMBER"

echo y | lpad reset

DESC="2 generations completion test." WC_ANALYZE_FAST=1 SINGLE_DAUGHTERS=1 N_GENS=2 MASS_DISTRIBUTION=0 \
	PARALLEL_PARCA=1 COMPRESS_OUTPUT=0 PLOTS=ACTIVE BUILD_CAUSALITY_NETWORK=1 RAISE_ON_TIME_LIMIT=1 python runscripts/fireworks/fw_queue.py

bash runscripts/jenkins/run-fireworks.sh

# Test some runscripts to make sure they do not fall out of step with other updates
# Remove these lines if these scripts are no longer useful and causing build failures
runscripts/debug/metabolism.py --validation -1  # This should come before parameter_search.py to use the correct default directory
runscripts/manual/parameter_search.py --method quick_example

rm -fr out/*
