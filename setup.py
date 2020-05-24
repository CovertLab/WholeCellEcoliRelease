from __future__ import absolute_import, division, print_function

from distutils.core import setup# , Extension
# from distutils.sysconfig import get_python_inc

from Cython.Build import cythonize

import numpy as np
import os

# The commented out code refers to files that have been removed; I've retained
# the code for future reference.

# pythonIncDir = get_python_inc()
# numpyIncDir = os.path.join(os.path.dirname(np.__file__), "core", "include", "numpy")

# polymerize_module = Extension(
# 	name = "wholecell.utils._polymerize",
# 	sources = ["wholecell/utils/polymerize.c"],
# 	include_dirs = [pythonIncDir, numpyIncDir],
# 	libraries = ["gsl", "gslcblas"],
# 	extra_compile_args = ["-fPIC"]
# 	)

# setup(name = "Polymerize",
# 	version = "0.0.1",
# 	description = "Polymerize module",
# 	ext_modules = [polymerize_module]
# 	)

build_sequences_module = cythonize(
	os.path.join("wholecell", "utils", "_build_sequences.pyx"),
	# annotate=True,
	)

setup(
	name = "Build sequences",
	ext_modules = build_sequences_module,
	include_dirs = [np.get_include()]
	)

complexation_module = cythonize(
	os.path.join("wholecell", "utils", "mc_complexation.pyx"),
	# annotate=True,
	)

setup(
	name = "Monte-carlo complexation",
	ext_modules = complexation_module,
	include_dirs = [np.get_include()]
	)

fast_polymerize_sums_module = cythonize(
	os.path.join("wholecell", "utils", "_fastsums.pyx"),
	#compiler_directives = {'linetrace': True},
	# annotate=True, # emit an html file with annotated C code
	)

setup(
	name = "Fast polymerize sums",
	ext_modules = fast_polymerize_sums_module,
	include_dirs = [np.get_include()]
	)
