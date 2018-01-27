"""
Test NumPy performance to ensure the configuration is good.

Running it this way prints all timing measurements:
	python -m wholecell.tests.utils.test_numpy_performance

Running it this way only prints timing measurements for failed tests,
e.g. those that exceed their @nose.tools.timed() thresholds:
	nosetests wholecell/tests/utils/test_numpy_performance.py
"""

import math
import resource
import time
import unittest

import numpy as np
import nose.plugins.attrib as noseAttrib
import nose.tools


def clock2():
	"""
	clock2() -> (ru_utime, ru_stime)

	Return the TOTAL (USER, SYSTEM CPU) time in floating point seconds
	since the start of the process, excluding any child processes. The
	data comes from resource.getrusage() so it avoids the wraparound
	problems in time.clock().

	Note: getrusage() also returns memory usage data.

	Code cribbed from IPython %time magic
	github.com/ipython/ipython - IPython/core/magics/execution.py
	"""
	return resource.getrusage(resource.RUSAGE_SELF)[:2]

def _format_time(timespan, precision=3):
	"""
	Formats the timespan in a human-readable form.

	Code cribbed from IPython %time magic
	github.com/ipython/ipython - IPython/core/magics/execution.py
	"""

	if timespan >= 60.0:
		# More than a minute. Format it in a human readable form.
		parts = [("d", 60*60*24), ("h", 60*60), ("min", 60), ("s", 1)]
		time = []
		leftover = timespan
		for suffix, length in parts:
			value = int(leftover / length)
			if value > 0:
				leftover = leftover % length
				time.append(u'%s%s' % (str(value), suffix))
			if leftover < 1:
				break
		return " ".join(time)

	units = [u"s", u"ms", u'us', "ns"]  # the save value
	scaling = [1, 1e3, 1e6, 1e9]

	if timespan > 0.0:
		order = min(-int(math.floor(math.log10(timespan)) // 3), 3)
	else:
		order = 3
	return u"%.*g %s" % (precision, timespan * scaling[order], units[order])


class Test_numpy_performance(unittest.TestCase):

	def setUp(self):
		self.M = np.random.random(size=(1000, 1000))

	def time_it(self, code_to_measure):
		"""
		Times the execution of code_to_measure().
		Code cribbed from IPython %time magic
		github.com/ipython/ipython - IPython/core/magics/execution.py
		"""
		wtime = time.time

		# time execution
		wall_st = wtime()
		st = clock2()
		code_to_measure()
		end = clock2()
		wall_end = wtime()

		wall_time = wall_end - wall_st
		cpu_user = end[0] - st[0]
		cpu_sys = end[1] - st[1]
		cpu_tot = cpu_user + cpu_sys
		print("%s CPU times: user %s, sys: %s, total: %s; Wall time: %s"
			% (self.id(), _format_time(cpu_user), _format_time(cpu_sys),
			   _format_time(cpu_tot), _format_time(wall_time)))

	@noseAttrib.attr('smalltest')
	@nose.tools.timed(0.050)
	def test_init(self):
		"""Time a NumPy matrix dot() operation."""
		M = self.M
		self.time_it(lambda: M.dot(M))


if __name__ == '__main__':
	unittest.main()
