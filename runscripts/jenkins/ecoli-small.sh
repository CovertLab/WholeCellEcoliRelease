set -e

module load wcEcoli/sherlock2
pyenv local wcEcoli-paper

make clean
make compile

PYTHONPATH=$PWD:$PYTHONPATH nosetests -a 'smalltest' --with-xunit --with-coverage --cover-package=wholecell --cover-xml
