"""
Test library performance (NumPy and libraries above and below it) to
discover configuration problems such as performance bugs in some versions
of pip packages or problems linking to native libraries. Precise timings
aren't needed.

Running it this way prints all timing measurements:
	python -m wholecell.tests.utils.test_numpy_performance

Running it this way prints timing measurements (and other printout) only
for failed tests, e.g. those that exceed their @nose.tools.timed()
thresholds:
	nosetests wholecell/tests/utils/test_numpy_performance.py

Running it this way runs the iterative test that isn't automatically
discovered as a test method:
	python -m unittest -v wholecell.tests.utils.test_numpy_performance.Test_numpy_performance.multitest_dot
"""

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

	FYI: getrusage() also returns memory usage data.

	Cf. IPython %time magic github.com/ipython/ipython -
	IPython/core/magics/execution.py
	"""
	usage = resource.getrusage(resource.RUSAGE_SELF)
	return usage.ru_utime, usage.ru_stime

def _format_time(timespan, precision=3):
	"""
	Formats the timespan in a human-readable form.

	Cf. IPython %time magic
	github.com/ipython/ipython - IPython/core/magics/execution.py
	"""

	if timespan >= 60.0:
		# More than a minute. Format it in a human readable form.
		PARTS = (("d", 60 * 60 * 24), ("h", 60 * 60), ("min", 60), ("s", 1))
		time_parts = []
		leftover = timespan
		for (suffix, length) in PARTS:
			value = int(leftover / length)
			if value > 0:
				leftover = leftover % length
				time_parts.append(u'%s%s' % (str(value), suffix))
			if leftover < 1:
				break
		return " ".join(time_parts)

	UNITS = ("s", "ms", "us", "ns")
	SCALING = (1, 1e3, 1e6, 1e9)
	K = 3  # orders of magnitude between SCALING factors

	if timespan > 0.0:
		order = min(-int(np.floor(np.log10(timespan)) // K), K)
	else:
		order = K
	return "%.*g %s" % (precision, timespan * SCALING[order], UNITS[order])


def time_it(code_to_measure, title='Measured'):
	"""
	Times the execution of code_to_measure().
	Cf. IPython %time magic
	github.com/ipython/ipython - IPython/core/magics/execution.py time()
	"""
	wall_clock = time.time

	# time execution
	wall_start = wall_clock()
	(start_user, start_sys) = clock2()
	code_to_measure()
	(end_user, end_sys) = clock2()
	wall_end = wall_clock()

	wall_time = wall_end - wall_start
	cpu_user = end_user - start_user
	cpu_sys = end_sys - start_sys
	cpu_total = cpu_user + cpu_sys

	print("%s CPU times: user %s, sys: %s, total: %s; Wall time: %s"
		  % (title, _format_time(cpu_user), _format_time(cpu_sys),
			 _format_time(cpu_total), _format_time(wall_time)))


class Test_numpy_performance(unittest.TestCase):
	"""Test one or more NumPy array ops to see that it's performing OK."""

	def time_this(self, code_to_measure):
		"""Times the execution of code_to_measure()."""
		time_it(code_to_measure, self.id())

	# On 2015 MacBook Pro this takes < 25 ms.
	# On Sherlock 1.0 with 1 CPU this takes ~250 ms.
	# On Sherlock 1.0 with 16 CPUs this takes 50 - 100 ms.
	# Allow time for test framework overhead + matrix construction.
	@noseAttrib.attr('smalltest')
	@nose.tools.timed(0.35)
	def test_dot(self):
		"""Time NumPy matrix dot()."""
		M = np.random.random(size=(1000, 1000))
		self.time_this(lambda : M.dot(M))

	def multitest_dot(self):
		"""Time NumPy matrix dot() many times."""
		for iteration in xrange(100):
			self.test_dot()


if __name__ == '__main__':
	unittest.main()
