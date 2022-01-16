set -e

source runscripts/jenkins/setup-environment.sh

pytest --cov=wholecell --cov-report xml --junitxml=unittests.xml
mypy
