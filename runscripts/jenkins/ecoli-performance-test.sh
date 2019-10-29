module load wcEcoli/sherlock2
pyenv local wcEcoli-paper

make clean
make compile

set -e

PYTHONPATH=$PWD:$PYTHONPATH nosetests -a 'performance' --with-xunit --with-coverage --cover-package=wholecell --cover-xml
