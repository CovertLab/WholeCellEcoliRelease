set -e

export PYTHONPATH=$PWD
module load wcEcoli/python3

### Edit this line to make this branch use another pyenv like wcEcoli3-staging
pyenv local wcEcoli3

make clean compile
