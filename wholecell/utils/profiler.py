'''
File for profiling tools
'''

from line_profiler import LineProfiler
from functools import wraps

def line_profile(func):
	'''
	Decorator for line profiling of a given function.
	Prints time spent on each line within a function every time the
	function is called.

	@profile decorator from LineProfiler requires running the script
	through kernprof and additional setup but this decorator can be
	used in front of any function in the model with the lines below.

	Usage:
		from wholecell.utils.profiler import line_profile
		@line_profile
		def function_to_profile():
			...
	'''

	@wraps(func)
	def wrapper(*args, **kwargs):
		lp = LineProfiler()
		try:
			return lp(func)(*args, **kwargs)
		finally:
			lp.print_stats()

	return wrapper
