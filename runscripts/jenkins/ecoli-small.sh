set -e

module load wcEcoli/python3
pyenv local wcEcoli3

make clean compile

PYTHONPATH=$PWD:$PYTHONPATH pytest --cov=wholecell --cov-report xml \
    --junitxml=unittests.xml

runscripts/debug/mypy.sh
