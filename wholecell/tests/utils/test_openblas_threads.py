#!/usr/bin/env python
"""
Any of these ways works to run this test. The first one is quiet if the test
passes. The others always show stdout:
	pytest wholecell/tests/utils/test_openblas_threads.py
	python -m wholecell.tests.utils.test_openblas_threads
	wholecell/tests/utils/test_openblas_threads.py
"""

import os
import sys
import unittest

import numpy as np

import wholecell.utils.filepath as fp
from wholecell.utils import parallelization


THIS_DIR_PATH = os.path.dirname(os.path.realpath(__file__))


def dot_product():
	"""Return a dot product of 2 data arrays. In some installations, the result
	varies with the number of OpenBLAS threads!
	"""
	diff = np.load(os.path.join(THIS_DIR_PATH, 'diff.npy'))
	mass = np.load(os.path.join(THIS_DIR_PATH, 'mass.npy'))
	return diff.dot(mass)


class Test_openblas_threads(unittest.TestCase):
	def test_openblas(self):
		"""Compare a dot product running with various numbers of OpenBLAS threads."""
		products = []
		thread_range = list(range(1, parallelization.cpus(advice='mac override') + 1)) + ['']
		print('{:>7} {:>26} {:>26}'.format('THREADS', 'DOT PRODUCT', 'DIFF FROM 1 THREAD'))

		for num_threads in thread_range:
			env = dict(os.environ, OPENBLAS_NUM_THREADS=str(num_threads))
			command = 'python -m wholecell.tests.utils.test_openblas_threads DOT'.split()
			dot = float(fp.run_cmd(command, env=env))

			products.append(dot)
			diff = dot - products[0]
			print('{:7} {:26.17g} {:26.17g}'.format(num_threads, dot, diff))

		assert products[0] == 0.016683805584112754

		# Check that multi-threaded results are within some tolerance.
		np_products = np.array(products)
		reference = np.full(len(products), products[0])
		np.testing.assert_allclose(np_products, reference, rtol=1.5e-14)


if __name__ == '__main__':
	if sys.argv[1:] == ['DOT']:
		print(dot_product())
	else:
		unittest.main()
