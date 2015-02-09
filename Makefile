.PHONY: compile, runSimulation, runSimulationJob, runAnalysisSingle, justKb, justSimulation, buildKb, fitKb_1, execModel_1, clean, clobber

compile:
	python2.7 setup.py build_ext --inplace
	rm -fr build

runSimulation: compile
	PYTHONPATH="${PWD}:${PYTHONPATH}" ./runscripts/runSimulation.sh

NSIMS ?= "4"

runSimulationJob: compile
	PYTHONPATH="${PWD}:${PYTHONPATH}" ./runscripts/queueSimulationAndAnalysis.sh $(NSIMS)

FIXTURES_KBDIR ?= "fixtures/kb"
FIXTURES_SIMDIR ?= "fixtures/sim"

buildKb: compile
	PYTHONPATH="${PWD}:${PYTHONPATH}" python2.7 runscripts/buildKb.py $(FIXTURES_KBDIR)

fitKb_1: compile
	PYTHONPATH="${PWD}:${PYTHONPATH}" python2.7 runscripts/fit.py 1 $(FIXTURES_KBDIR) $(FIXTURES_SIMDIR)

execModel_1: compile
	PYTHONPATH="${PWD}:${PYTHONPATH}" python2.7 runscripts/execModel.py 1 $(FIXTURES_KBDIR) $(FIXTURES_SIMDIR)

execAnalysis_1: compile
	runscripts/execAnalysis.sh 1 $(FIXTURES_KBDIR) $(FIXTURES_SIMDIR)

execAnalysis_2: compile
	runscripts/execAnalysis.sh 2 $(FIXTURES_KBDIR) $(FIXTURES_SIMDIR)

fitKb_2: compile
	PYTHONPATH="${PWD}:${PYTHONPATH}" python2.7 runscripts/fit.py 2 $(FIXTURES_KBDIR) $(FIXTURES_SIMDIR)

execModel_2: compile
	PYTHONPATH="${PWD}:${PYTHONPATH}" python2.7 runscripts/execModel.py 2 $(FIXTURES_KBDIR) $(FIXTURES_SIMDIR)

justKb: buildKb fitKb_1 # execModel_1 execAnalysis_1 fitKb_2

justSimulation: execModel_2

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
