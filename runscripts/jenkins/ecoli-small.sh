set -e

module load wcEcoli/sherlock2
pyenv local wcEcoli2

make clean
make compile

PYTHONPATH=$PWD:$PYTHONPATH pytest --cov=wholecell --cov-report xml \
    --junitxml=unittests.xml
