module load wcEcoli/python3
pyenv local wcEcoli3

make clean
make compile

set -e

# Running it this way prints all timing measurements:
python -m wholecell.tests.utils.test_library_performance
