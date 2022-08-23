set -e

export PYTHONPATH=$PWD
module load wcEcoli/python3

if [ -d "${PYENV_ROOT}" ]; then
    export PATH="${PYENV_ROOT}/bin:${PATH}"
    eval "$(pyenv init -)"
    eval "$(pyenv virtualenv-init -)"
fi

### Edit this line to make this branch use another pyenv like wcEcoli3-staging
pyenv local wcEcoli3
pyenv activate

make clean compile
