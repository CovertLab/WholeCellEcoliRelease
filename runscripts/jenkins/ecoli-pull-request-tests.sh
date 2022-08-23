set -e

### ---------------------------------------------------------------------------
### Edit the pyenv in setup-environment.sh to make the PR build temporarily use
### another pyenv for testing (eg. wcEcoli3-staging). Revert the change before
### merging the PR into master to prevent changing it for other Jenkins builds.
### ---------------------------------------------------------------------------
source runscripts/jenkins/setup-environment.sh

echo pyenv version: $(pyenv version)

pytest --cov=wholecell --cov-report xml --junitxml=unittests.xml
