.PHONY: compile, runSimulation, runSimulationJob, runAnalysisSingle, justKb, justSimulation, buildKb, clean, clobber

compile:
	python2.7 setup.py build_ext --inplace
	rm -fr build

runSimulation: compile
	PYTHONPATH="${PWD}:${PYTHONPATH}" ./runscripts/runSimulation.sh

runSimulationJob: compile
	PYTHONPATH="${PWD}:${PYTHONPATH}" ./runscripts/queueSimulationAndAnalysis.sh 4

# TODO: Get rid of this target?
runAnalysisSingle:
	./runscripts/runAnalysisSingle.sh out/simOut out/plotOut wholecell/analysis/single/

justKb: compile
	PYTHONPATH="${PWD}:${PYTHONPATH}" python2.7 runscripts/createKbs.py

justSimulation: compile
	PYTHONPATH="${PWD}:${PYTHONPATH}" python2.7 runscripts/justSimulation.py

buildKb: compile
	PYTHONPATH="${PWD}:${PYTHONPATH}" python2.7 runscripts/buildKb.py

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
