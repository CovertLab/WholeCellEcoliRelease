.PHONY: compile, clean

compile:
	python2.7 setup.py build_ext --inplace
	rm -fr build

clean:
	rm -fr fixtures
	(cd reconstruction/ecoli/dataclasses/process && rm -f equilibrium_odes.py two_component_system_odes*.py)
	find . -not \( -path ./out -prune \) -not \( -path ./.git -prune \) -name "*.cPickle" -exec rm -fr {} \;
	find . -name "*.pyc" -exec rm -rf {} \;
	find . -name "*.o" -exec rm -fr {} \;
	find . -name "*.so" -exec rm -fr {} \;
	rm -fr build
	rm -fr launcher_20* block_20*
