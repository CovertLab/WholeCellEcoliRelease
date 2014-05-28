.PHONY: all, runSimulation, runSimulationJob, runAnalysisSingle, clean, clobber

all: wholecell/utils/_polymerize.so

wholecell/utils/_polymerize.so: wholecell/utils/polymerize.c wholecell/utils/polymerize.h
	python2.7 setup.py build_ext --inplace
	rm -fr build

runSimulation: all
	./runscripts/runSimulation.sh

runSimulationJob: all
	./runscripts/queueSimulationAndAnalysis.sh 4

runAnalysisSingle:
	./runscripts/runAnalysisSingle.sh out/simOut out/plotOut wholecell/analysis/single/

clean:
	find . -name "*.cPickle" -exec rm -fr {} \;
	find . -name "*.pyc" -exec rm -rf {} \;
	find . -name "*.o" -exec rm -fr {} \;
	find . -name "*.so" -exec rm -fr {} \;
	rm -fr build

clobber:
	rm -fr out/simOut
	rm -fr out/plotOut
	rm -fr runAnalysis*
	rm -fr runSimulation*
	rm -fr simShellLog*
	rm -fr analysisSingle*

### Targets for tests are below

.PHONY: testFixtures

testFixtures:
	python2.7 fixtures/generateTestFixtures.py