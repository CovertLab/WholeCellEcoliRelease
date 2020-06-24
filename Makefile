.PHONY: compile, clean, recompile

compile:
	python setup.py build_ext --inplace
	rm -fr build

# The reconstruction/ecoli/dataclasses/process/*.py files were generated by
#  write_ode_file.py in Parca code.
# Fireworks writes launcher_20* and block_20*.
clean:
	rm -fr fixtures
	(cd reconstruction/ecoli/dataclasses/process && rm -f equilibrium_odes.py two_component_system_odes*.py)
	find . -name "*.pyc" -exec rm -rf {} \;
	find . -name "*.o" -exec rm -fr {} \;
	find . -name "*.so" -exec rm -fr {} \;
	rm -fr build
	rm -fr launcher_20* block_20*

# Delete just the *.so libraries then (re)compile them.
# This is useful when switching to a different Python virtualenv.
recompile:
	find . -name "*.so" -delete
	make compile
