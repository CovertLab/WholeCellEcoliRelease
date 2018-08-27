module load wcEcoli/sherlock2
pyenv local wcEcoli2

make clean
make compile

PYTHONPATH=$PWD:$PYTHONPATH nosetests -a 'performance' --with-xunit --with-coverage --cover-package=wholecell --cover-xml
