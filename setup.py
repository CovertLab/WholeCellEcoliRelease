#from setuptools import setup
from distutils.core import setup

from Cython.Build import cythonize

import numpy as np
import os


build_sequences_module = cythonize(
	os.path.join("wholecell", "utils", "_build_sequences.pyx"),
	# annotate=True,
	)

setup(
	name = "Build sequences",
	ext_modules = build_sequences_module,
	include_dirs = [np.get_include()],
	# zip_safe = False  # see Cython docs, but it seems to need setuptools
	)

complexation_module = cythonize(
	os.path.join("wholecell", "utils", "mc_complexation.pyx"),
	# annotate=True,
	)

setup(
	name = "Monte-carlo complexation",
	ext_modules = complexation_module,
	include_dirs = [np.get_include()],
	# zip_safe = False  # see Cython docs, but it seems to need setuptools
	)

fast_polymerize_sums_module = cythonize(
	os.path.join("wholecell", "utils", "_fastsums.pyx"),
	#compiler_directives = {'linetrace': True},
	# annotate=True, # emit an html file with annotated C code
	)

setup(
	name = "Fast polymerize sums",
	ext_modules = fast_polymerize_sums_module,
	include_dirs = [np.get_include()],
	# zip_safe = False  # see Cython docs, but it seems to need setuptools
	)
