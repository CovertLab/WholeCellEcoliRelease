set -e

module load wcEcoli/sherlock2
pyenv local wcEcoli2

make clean
make compile

PYTHONPATH=$PWD:$PYTHONPATH nosetests -a 'smalltest' --with-xunit --with-coverage --cover-package=wholecell --cover-xml
