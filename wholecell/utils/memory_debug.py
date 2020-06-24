"""
Memory leak debugging utility. This should detect uncollectable Python objects
but not leaked C/C++ nodes.

If the caller defaults `enabled`, then the OS environment variable `DEBUG_GC`
will enable/disable memory leak detection.

Example:
	with detect_leaks():
		make_a_matplotlib_chart()

`detect_leaks()` is a context manager that runs a full GC and restores GC debug
flags when done. It should work fine to nest `with detect_leaks(False)` for a
block of code inside `with detect_leaks(True)`.
"""

from __future__ import absolute_import, division, print_function

from contextlib import contextmanager
import gc
import os
import sys
import traceback


# See https://pymotw.com/2/gc/ for more info on using the gc interface.
# DEBUG_UNCOLLECTABLE causes the collector to report on objects it can't
# collect. You need to combine it with DEBUG_OBJECTS to print info about
# objects and DEBUG_INSTANCES to print info about instances of old-style
# classes (not derived from object).
#
# See https://docs.python.org/3/library/gc.html for Python 3.4+. Following
# PEP 442, objects with a __del__() method don't end up in gc.garbage anymore.
DEBUG_OBJECTS = getattr(gc, 'DEBUG_OBJECTS', 0)  # defined in python 2
DEBUG_INSTANCES = getattr(gc, 'DEBUG_INSTANCES', 0)  # defined in python 2
TRACE_UNCOLLECTABLES = gc.DEBUG_UNCOLLECTABLE | DEBUG_OBJECTS | DEBUG_INSTANCES
TRACE_NONE = 0


@contextmanager
def detect_leaks(enabled=None):
	"""A context manager that optionally detects Python object leaks in the
	`with` statement body.

	Set `enabled` to True to enable leak detection, False to disable leak
	detection, or default to let the 'DEBUG_GC' environment variable
	(int, "0") enable leak detection if non-zero.

	Leak detection has some overhead including running a full collection and
	printing a list of uncollectable objects.

	Per https://docs.python.org/2/library/gc.html, "Objects that have
	__del__() methods and are part of a reference cycle cause the entire
	reference cycle to be uncollectable, including objects not necessarily in
	the cycle but reachable only from it. Python doesn't collect such cycles
	automatically because, in general, it isn't possible for Python to guess
	a safe order in which to run the __del__() methods."
	"""
	if enabled is None:
		enabled = bool(int(os.environ.get('DEBUG_GC', '0')))

	saved_debug_flags = gc.get_debug()
	gc.set_debug(TRACE_UNCOLLECTABLES if enabled else TRACE_NONE)

	yield  # yield to the `with` statement body

	if enabled:
		if gc.collect():  # prints lines like "gc: uncollectable <CleanupGraph 0x10045f810>"
			# Print some stack trace to show where the uncollectables were found
			print('gc uncollectables detected in:', file=sys.stderr)
			traceback.print_stack(limit=3)
			# Examine the gc.garbage list here?

	gc.set_debug(saved_debug_flags)
