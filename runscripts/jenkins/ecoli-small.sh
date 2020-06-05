set -e

module load wcEcoli/sherlock2

WCECOLI_PYENV=wcEcoli2
pyenv local ${WCECOLI_PYENV}

make clean compile

PYTHONPATH=$PWD:$PYTHONPATH pytest --cov=wholecell --cov-report xml \
    --junitxml=unittests.xml

(export PYENV_VERSION="mypy:${WCECOLI_PYENV}"; echo ---Running mypy---; mypy --py2)
