"""
Abstract base class for a Firetask that runs a category of analysis plots.

If the `DEBUG_GC` environment variable is true, enable memory leak detection.
"""

from __future__ import absolute_import, division, print_function

import abc
from collections import OrderedDict
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

		TAGS = the dictionary that maps tag names to lists of analysis plot
			module filenames (in this category) to run.

			Use all uppercase for tag names so they don't conflict with module
			filenames.

			The tag 'CORE' lists the plots to run by default when the argument
			list is empty.

			The tag 'ACTIVE' lists all active plots in this category. The
			nightly build should run 'ACTIVE'.

		MODULE_PATH = the module pathname for loading this subclass's analysis
			plots.

	Optional params include plots_to_run, output_filename_prefix, cpus.
	"""

	@abc.abstractmethod
	def plotter_args(self, module_filename):
		"""(Abstract) Return a tuple of arguments to pass to the analysis plot
		class' `main()` method.
		"""
		raise NotImplementedError

	def expand_plot_names(self, plot_names, name_dict):
		'''Recursively expand TAGS and plot names that lack the '.py' suffix,
		adding the names to name_dict. name_dict is an OrderedDict doing the
		job of an OrderedSet class.
		'''
		for name in plot_names:
			if name in self.TAGS:
				self.expand_plot_names(self.TAGS[name], name_dict)
			elif name.endswith('.py'):
				name_dict[name] = True
			else:
				name_dict[name + '.py'] = True

	def list_plot_files(self, plot_names):
		'''List the plot module files (within self.MODULE_PATH) named by the
		given list of plot names, doing these transformations:

			* Default to the 'CORE' tag if plot_names is empty.
			* Expand all TAGS as defined by self.TAGS.
			* Append '.py' to filenames as needed.
			* Deduplicate entries but preserve the order.
		'''
		if not plot_names:
			plot_names = ['CORE']
		name_dict = OrderedDict()
		self.expand_plot_names(plot_names, name_dict)
		return name_dict.keys()

	def run_task(self, fw_spec):
		startTime = time.time()
		print("\n{}: --- Starting {} ---".format(
			time.ctime(startTime), type(self).__name__))

		plot_names = self.get("plots_to_run", [])
		fileList = self.list_plot_files(plot_names)

		self['output_filename_prefix'] = self.get('output_filename_prefix', '')

		cpus = min(self.get("cpus", 1), parallelization.cpus())
		pool = None
		results = {}

		if cpus > 1:
			pool = mp.Pool(processes=cpus)

		exceptionFileList = []
		for f in fileList:
			try:
				mod = importlib.import_module(self.MODULE_PATH + '.' + f[:-3])
			except ImportError:
				traceback.print_exc()
				exceptionFileList.append(f)
				continue

			args = self.plotter_args(f)

			if pool:
				results[f] = pool.apply_async(run_plot, args=(mod.Plot, args, f))
			else:
				print("{}: Running {}".format(time.ctime(), f))
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
			print("Completed analysis in {} with an exception in:".format(duration))
			for file in exceptionFileList:
				print("\t{}".format(file))
			raise Exception("Error in analysis")
		else:
			print("Completed analysis in {}".format(duration))


def run_plot(plot_class, args, name):
	"""Run the given plot class in a Pool worker.
	Since this Firetask is running multiple plot classes in parallel, ask them
	to use just 1 CPU core each.
	"""
	try:
		print("{}: Running {}".format(time.ctime(), name))
		plot_class.main(*args, cpus=1)
	except KeyboardInterrupt:
		sys.exit(1)
	except Exception as e:
		traceback.print_exc()
		raise Exception(e)  # TODO: Return e so the caller can print it?
