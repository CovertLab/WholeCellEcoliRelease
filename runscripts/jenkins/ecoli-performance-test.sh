source runscripts/jenkins/setup-environment.sh

set -e

# Running it this way prints all timing measurements:
python -m wholecell.tests.utils.test_library_performance
