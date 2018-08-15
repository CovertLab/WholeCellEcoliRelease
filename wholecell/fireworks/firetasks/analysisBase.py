"""
Abstract base class for a Firetask that runs a category of analysis plots.

If the `DEBUG_GC` environment variable is true, enable memory leak detection.
"""

from __future__ import absolute_import
from __future__ import division

import abc
import importlib
import multiprocessing as mp
import sys
import time
import traceback

from fireworks import FiretaskBase

from wholecell.utils import parallelization


class AnalysisBase(FiretaskBase):
	"""Base class for analysis plot Firetasks.

	Each subclass should set the usual Firetask class variables _fw_name,
	required_params, and optional_params; also

		ACTIVE_MODULES = the list of active module names for this category to
			use if the caller doesn't provide a specific list.
		MODULE_PATH = the module pathname to load for this subclass of analysis
			plots.

	Optional params include plots_to_run, output_filename_prefix, cpus.
	"""

	@abc.abstractmethod
	def plotter_args(self, module_filename):
		"""(Abstract) Return a tuple of arguments to pass to the analysis plot
		class' `main()` method.
		"""
		raise NotImplementedError

	def run_task(self, fw_spec):
		startTime = time.time()
		print "\n{}: --- Starting {} ---".format(
			time.ctime(startTime), type(self).__name__)

		fileList = self.get("plots_to_run", [])
		if not fileList:
			fileList = self.ACTIVE_MODULES

		self['output_filename_prefix'] = self.get('output_filename_prefix', '')

		cpus = min(self.get("cpus", 1), parallelization.cpus())
		pool = None
		results = {}

		if cpus > 1:
			pool = mp.Pool(processes=cpus)

		exceptionFileList = []
		for f in fileList:
			mod = importlib.import_module(self.MODULE_PATH + '.' + f[:-3])
			args = self.plotter_args(f)

			if pool:
				results[f] = pool.apply_async(run_plot, args=(mod.Plot, args, f))
			else:
				print "{}: Running {}".format(time.ctime(), f)
				try:
					mod.Plot.main(*args)
				except Exception:
					traceback.print_exc()
					exceptionFileList.append(f)

		if pool:
			pool.close()
			pool.join()
			for f, result in results.items():
				if not result.successful():
					exceptionFileList.append(f)

		timeTotal = time.time() - startTime

		duration = time.strftime("%H:%M:%S", time.gmtime(timeTotal))
		if exceptionFileList:
			print "Completed analysis in {} with an exception in:".format(duration)
			for file in exceptionFileList:
				print "\t{}".format(file)
			raise Exception("Error in analysis")
		else:
			print "Completed analysis in {}".format(duration)


def run_plot(plot_class, args, name):
	"""Run the given plot class in a Pool worker.
	Since this Firetask is running multiple plot classes in parallel, ask them
	to use just 1 CPU core each.
	"""
	try:
		print "{}: Running {}".format(time.ctime(), name)
		plot_class.main(*args, cpus=1)
	except KeyboardInterrupt:
		sys.exit(1)
	except Exception as e:
		traceback.print_exc()
		raise Exception(e)  # TODO: Return e so the caller can print it?
