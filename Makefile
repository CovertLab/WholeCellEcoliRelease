.PHONY: clean, all, runSimulation

all: wholecell/utils/_polymerize.so

wholecell/utils/_polymerize.so: wholecell/utils/polymerize.c wholecell/utils/polymerize.h
	python setup.py build_ext --inplace
	rm -fr build

runSimulation: all
	python runscripts/runSimulation.py

clean:
	find . -name "*.cPickle" -exec rm -fr {} \;
	find . -name "*.pyc" -exec rm -rf {} \;
	find . -name "*.o" -exec rm -fr {} \;
	find . -name "*.so" -exec rm -fr {} \;
	rm -fr build