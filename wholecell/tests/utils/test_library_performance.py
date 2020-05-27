"""
Test library performance (NumPy and libraries above and below it) to
discover configuration problems such as performance bugs in some versions
of pip packages or problems linking to native libraries. Precise timings
aren't needed.

Running it this way prints all timing measurements:
	python -m wholecell.tests.utils.test_library_performance

Running it these ways prints timing measurements (and other printout) only
for failed tests:
	pytest wholecell/tests/utils/test_library_performance.py

Running it this way runs the iterative test that isn't automatically
discovered as a test method:
	python -m unittest -v wholecell.tests.utils.test_library_performance.Test_library_performance.multitest_dot
"""

from __future__ import absolute_import, division, print_function

import resource
import time
import unittest

import numpy as np
import scipy.integrate


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
				time_parts.append("%s%s" % (str(value), suffix))
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
	Time code_to_measure(). Print timings. Return the elapsed wall time.
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

	print("\n%s CPU time: user %s + sys %s = %s; Wall: %s"
		  % (title, _format_time(cpu_user), _format_time(cpu_sys),
			 _format_time(cpu_total), _format_time(wall_time)))
	return wall_time


# Originally generated by the parameter calculator (parca) into
# reconstruction/ecoli/dataclasses/process/two_component_system_odes_parca.py
#
# TODO(jerry): Use smaller matrices for this test.
#
# TODO(jerry): In the original parca-generated code, remove the `**1.0`
# terms. In both, factor out common subexpressions and constants?
#
# TODO(jerry): Try JIT-compiling this with Numba.
def derivatives(y, t):
	return np.array([
		[-100000000.0*y[0]*y[6] + 500.0*y[3]*y[4]], [0], [0],
		[100000000.0*y[0]*y[6] - 500.0*y[3]*y[4]], [0],
		[100000000.0*y[0]*y[6] - 0.01*y[5]*y[8]],
		[-100000000.0*y[0]*y[6] + 0.01*y[5]*y[8]], [0], [0],
		[170000.0*y[10]*y[4] - 100000000.0*y[12]*y[9]],
		[-170000.0*y[10]*y[4] + 100000000.0*y[12]*y[9]],
		[-0.01*y[11]*y[8] + 100000000.0*y[12]*y[13] + 100000000.0*y[12]*y[9]],
		[0.01*y[11]*y[8] - 100000000.0*y[12]*y[13] - 100000000.0*y[12]*y[9]],
		[-100000000.0*y[12]*y[13] + 0.0001*y[14]*y[4]],
		[100000000.0*y[12]*y[13] - 0.0001*y[14]*y[4]],
		[-100000000.0*y[15]*y[18] + 170000.0*y[16]*y[4]],
		[100000000.0*y[15]*y[18] - 170000.0*y[16]*y[4]],
		[100000000.0*y[15]*y[18] - 0.01*y[17]*y[8] + 100000000.0*y[18]*y[19]],
		[-100000000.0*y[15]*y[18] + 0.01*y[17]*y[8] - 100000000.0*y[18]*y[19]],
		[-100000000.0*y[18]*y[19] + 0.0001*y[20]*y[4]],
		[100000000.0*y[18]*y[19] - 0.0001*y[20]*y[4]],
		[-100000000.0*y[21]*y[24] + 170000.0*y[22]*y[4]],
		[100000000.0*y[21]*y[24] - 170000.0*y[22]*y[4]],
		[100000000.0*y[21]*y[24] - 0.01*y[23]*y[8] + 100000000.0*y[24]*y[25]],
		[-100000000.0*y[21]*y[24] + 0.01*y[23]*y[8] - 100000000.0*y[24]*y[25]],
		[-100000000.0*y[24]*y[25] + 0.0001*y[26]*y[4]],
		[100000000.0*y[24]*y[25] - 0.0001*y[26]*y[4]],
		[-100000000.0*y[27]*y[30] + 170000.0*y[28]*y[4]],
		[100000000.0*y[27]*y[30] - 170000.0*y[28]*y[4]],
		[100000000.0*y[27]*y[30] - 0.01*y[29]*y[8] + 100000000.0*y[30]*y[31]],
		[-100000000.0*y[27]*y[30] + 0.01*y[29]*y[8] - 100000000.0*y[30]*y[31]],
		[-100000000.0*y[30]*y[31] + 0.0001*y[32]*y[4]],
		[100000000.0*y[30]*y[31] - 0.0001*y[32]*y[4]],
		[-100000000.0*y[33]*y[36] + 500.0*y[34]*y[4]],
		[100000000.0*y[33]*y[36] - 500.0*y[34]*y[4]],
		[100000000.0*y[33]*y[36] - 0.01*y[35]*y[8]],
		[-100000000.0*y[33]*y[36] + 0.01*y[35]*y[8]],
		[-100000000.0*y[37]*y[40] + 500.0*y[38]*y[4]],
		[100000000.0*y[37]*y[40] - 500.0*y[38]*y[4]],
		[100000000.0*y[37]*y[40] - 0.01*y[39]*y[8]],
		[-100000000.0*y[37]*y[40] + 0.01*y[39]*y[8]]]).reshape(-1)

# Ditto.
def derivativesJacobian(y, t):
	return np.array([
		[-100000000.0*y[6], 0, 0, 500.0*y[4], 500.0*y[3], 0, -100000000.0*y[0], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[100000000.0*y[6], 0, 0, -500.0*y[4], -500.0*y[3], 0, 100000000.0*y[0], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[100000000.0*y[6], 0, 0, 0, 0, -0.01*y[8], 100000000.0*y[0], 0, -0.01*y[5], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[-100000000.0*y[6], 0, 0, 0, 0, 0.01*y[8], -100000000.0*y[0], 0, 0.01*y[5], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 170000.0*y[10], 0, 0, 0, 0, -100000000.0*y[12], 170000.0*y[4], 0, -100000000.0*y[9], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, -170000.0*y[10], 0, 0, 0, 0, 100000000.0*y[12], -170000.0*y[4], 0, 100000000.0*y[9], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0, 0, 0, 0, -0.01*y[11], 100000000.0*y[12], 0, -0.01*y[8], 100000000.0*y[13] + 100000000.0*y[9], 100000000.0*y[12], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0, 0, 0, 0, 0.01*y[11], -100000000.0*y[12], 0, 0.01*y[8], -100000000.0*y[13] - 100000000.0*y[9], -100000000.0*y[12], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0.0001*y[14], 0, 0, 0, 0, 0, 0, 0, -100000000.0*y[13], -100000000.0*y[12], 0.0001*y[4], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, -0.0001*y[14], 0, 0, 0, 0, 0, 0, 0, 100000000.0*y[13], 100000000.0*y[12], -0.0001*y[4], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 170000.0*y[16], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -100000000.0*y[18], 170000.0*y[4], 0, -100000000.0*y[15], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, -170000.0*y[16], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100000000.0*y[18], -170000.0*y[4], 0, 100000000.0*y[15], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0, 0, 0, 0, -0.01*y[17], 0, 0, 0, 0, 0, 0, 100000000.0*y[18], 0, -0.01*y[8], 100000000.0*y[15] + 100000000.0*y[19], 100000000.0*y[18], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0, 0, 0, 0, 0.01*y[17], 0, 0, 0, 0, 0, 0, -100000000.0*y[18], 0, 0.01*y[8], -100000000.0*y[15] - 100000000.0*y[19], -100000000.0*y[18], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0.0001*y[20], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -100000000.0*y[19], -100000000.0*y[18], 0.0001*y[4], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, -0.0001*y[20], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100000000.0*y[19], 100000000.0*y[18], -0.0001*y[4], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 170000.0*y[22], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -100000000.0*y[24], 170000.0*y[4], 0, -100000000.0*y[21], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, -170000.0*y[22], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100000000.0*y[24], -170000.0*y[4], 0, 100000000.0*y[21], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0, 0, 0, 0, -0.01*y[23], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100000000.0*y[24], 0, -0.01*y[8], 100000000.0*y[21] + 100000000.0*y[25], 100000000.0*y[24], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0, 0, 0, 0, 0.01*y[23], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -100000000.0*y[24], 0, 0.01*y[8], -100000000.0*y[21] - 100000000.0*y[25], -100000000.0*y[24], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0.0001*y[26], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -100000000.0*y[25], -100000000.0*y[24], 0.0001*y[4], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, -0.0001*y[26], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100000000.0*y[25], 100000000.0*y[24], -0.0001*y[4], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 170000.0*y[28], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -100000000.0*y[30], 170000.0*y[4], 0, -100000000.0*y[27], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, -170000.0*y[28], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100000000.0*y[30], -170000.0*y[4], 0, 100000000.0*y[27], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0, 0, 0, 0, -0.01*y[29], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100000000.0*y[30], 0, -0.01*y[8], 100000000.0*y[27] + 100000000.0*y[31], 100000000.0*y[30], 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0, 0, 0, 0, 0.01*y[29], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -100000000.0*y[30], 0, 0.01*y[8], -100000000.0*y[27] - 100000000.0*y[31], -100000000.0*y[30], 0, 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 0.0001*y[32], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -100000000.0*y[31], -100000000.0*y[30], 0.0001*y[4], 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, -0.0001*y[32], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100000000.0*y[31], 100000000.0*y[30], -0.0001*y[4], 0, 0, 0, 0, 0, 0, 0, 0],
		[0, 0, 0, 0, 500.0*y[34], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -100000000.0*y[36], 500.0*y[4], 0, -100000000.0*y[33], 0, 0, 0, 0],
		[0, 0, 0, 0, -500.0*y[34], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100000000.0*y[36], -500.0*y[4], 0, 100000000.0*y[33], 0, 0, 0, 0],
		[0, 0, 0, 0, 0, 0, 0, 0, -0.01*y[35], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100000000.0*y[36], 0, -0.01*y[8], 100000000.0*y[33], 0, 0, 0, 0],
		[0, 0, 0, 0, 0, 0, 0, 0, 0.01*y[35], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -100000000.0*y[36], 0, 0.01*y[8], -100000000.0*y[33], 0, 0, 0, 0],
		[0, 0, 0, 0, 500.0*y[38], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -100000000.0*y[40], 500.0*y[4], 0, -100000000.0*y[37]],
		[0, 0, 0, 0, -500.0*y[38], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100000000.0*y[40], -500.0*y[4], 0, 100000000.0*y[37]],
		[0, 0, 0, 0, 0, 0, 0, 0, -0.01*y[39], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100000000.0*y[40], 0, -0.01*y[8], 100000000.0*y[37]],
		[0, 0, 0, 0, 0, 0, 0, 0, 0.01*y[39], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -100000000.0*y[40], 0, 0.01*y[8], -100000000.0*y[37]]])


class Test_library_performance(unittest.TestCase):
	"""Test some library operations to see that they're performing OK."""

	def time_this(self, code_to_measure, limit=0.5):
		"""Time the execution of code_to_measure(), enforcing a limit."""
		test_method_name = self.id().rpartition('.')[-1]
		elapsed = time_it(code_to_measure, test_method_name)
		message = "{} sec elapsed with {} sec limit".format(elapsed, limit)
		self.assertLessEqual(elapsed, limit, message)

	# On 2015 MacBook Pro this takes < 25 ms.
	# Sherlock 1.0 performance varies widely with number of CPUs,
	# OPENBLAS_NUM_THREADS=... value, compute node, and BLAS library.
	# Allow time for test framework overhead + matrix construction.
	def test_dot(self):
		"""Time NumPy float64 x float64 matrix dot()."""
		M = np.random.random(size=(1000, 1000))
		self.time_this(lambda: M.dot(M), 0.3)

	def multitest_dot(self):
		"""Time NumPy matrix dot() many times."""
		for iteration in xrange(100):
			self.test_dot()

	def test_int_dot_int(self):
		"""Time NumPy int64 x int64 matrix dot()."""
		N = np.random.randint(0, 10, size=(1000, 1000))
		self.time_this(lambda: N.dot(N), 5.0)  # SLOW!

	def test_int_dot_floated_int(self):
		"""
		Time converting an int64 matrix to float64 then int64 x float64
		matrix dot(). This is 30x - 90x faster than int64 x int64 because
		(1) modern CPUs have high-throughput floating point hardware,
		(2) BLAS has no integer type, and
		(3) the libraries don't parallelize integer matrix multiply.
		"""
		N = np.random.randint(0, 10, size=(1000, 1000))
		self.time_this(lambda: N.dot(N * 1.0), 0.3)

	@unittest.skip('pretty much the same as test_int_dot_floated_int()')
	def test_floated_int_dot_int(self):
		"""Time NumPy integer x float(integer) matrix dot()."""
		N = np.random.randint(0, 10, size=(1000, 1000))
		self.time_this(lambda: (N * 1.0).dot(N), 0.3)

	def test_int_dot_float(self):
		"""Time NumPy integer x float matrix dot()."""
		N = np.random.randint(0, 10, size=(1000, 1000))
		M = np.random.random(size=(1000, 1000))
		self.time_this(lambda: N.dot(M), 0.3)

	@unittest.skip('pretty much the same as test_int_dot_float()')
	def test_float_dot_int(self):
		"""Time NumPy float x integer matrix dot()."""
		M = np.random.random(size=(1000, 1000))
		N = np.random.randint(0, 10, size=(1000, 1000))
		self.time_this(lambda: M.dot(N), 0.3)

	def test_int_to_float32_dot_and_back(self):
		"""Time NumPy integer matrix converted to float32, dot(), and
		back. This can be twice as fast as float (float64) math.
		"""
		N = np.random.randint(0, 10000, size=(1000, 1000))
		M = np.random.random(size=(1000, 1000))
		self.time_this(lambda: N.astype(np.float32)
					   .dot(M.astype(np.float32)).astype(np.float32),
					   0.6)

	# Allow time for test framework overhead + matrix construction.
	def test_odeint(self):
		"""Time scipy.integrate.odeint()."""
		y0 = np.random.random(41)

		def odeint():
			y = scipy.integrate.odeint(
				derivatives, y0, t=[0, 1e6], Dfun=derivativesJacobian,
				mxstep=10000)
		self.time_this(odeint, 0.4)


if __name__ == '__main__':
	unittest.main()
