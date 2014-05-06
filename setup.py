from distutils.core import setup, Extension
from distutils.sysconfig import get_python_inc

import numpy as np
import os

pythonIncDir = get_python_inc()
numpyIncDir = os.path.join(os.path.dirname(np.__file__), "core", "include", "numpy")

polymerize_module = Extension(
	name = "wholecell.utils._polymerize",
	sources = ["wholecell/utils/polymerize.c"],
	include_dirs = [pythonIncDir, numpyIncDir],
	libraries = ["gsl", "gslcblas"],
	extra_compile_args = ["-fPIC"]
	)

setup(name = "Polymerize",
	version = "0.0.1",
	description = "Polymerize module",
	ext_modules = [polymerize_module]
	)
