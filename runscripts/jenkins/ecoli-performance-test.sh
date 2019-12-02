module load wcEcoli/sherlock2
pyenv local wcEcoli2

make clean
make compile

set -e

# Running it this way prints all timing measurements:
python -m wholecell.tests.utils.test_library_performance
