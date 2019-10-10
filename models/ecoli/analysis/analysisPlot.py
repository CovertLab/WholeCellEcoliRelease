"""
Common code for analysis plots. The abstract base class AnalysisPlot defines a
plot() method for scripts to call.

TODO: See Issue #161 Matplotlib backend enforcement. If updating to matplotlib
2.2.2 doesn't fix it, this class can implement a workaround but must do so
before any code imports pyplot, which means every subclass must import this
file before importing pyplot.

TODO: Reliably load the wcEcoli/matplotlibrc file even if the working directory
is wrong (see Issue #132). Setting the working directory could work if done
before matplotlib loads it, and that should also fix the backend (#161),
otherwise loading matplotlibrc here just ensures loading of other resources.

TODO: Enable future warnings, esp. for matplotlib.

TODO: Move the run_plot() args to instance variables?

TODO: Other shared code to simplify the subclasses, e.g. make plotOutDir,
check that `os.path.isdir(simOutDir)`, instantiate an AnalysisPaths (except for
SingleAnalysisPlot subclasses), etc.
"""

from __future__ import absolute_import, division, print_function

import abc
import matplotlib as mp
from matplotlib import pyplot as plt
from wholecell.utils import memory_debug, parallelization
from wholecell.utils import filepath as fp


class AnalysisPlot(object):
	"""Abstract Base Class for analysis plots.

	Each analysis class must override do_plot().

	Call main() to run an analysis plot for a Firetask.

	Use the environment variable 'DEBUG_GC' to enable memory leak debugging.

	Inputs:
		cpus: allotted number of CPU cores; default (0) => all available cores
	"""
	__metaclass__ = abc.ABCMeta

	def __init__(self, cpus=0):
		mp.rc_file(fp.MATPLOTLIBRC_FILE)
		self.cpus = parallelization.cpus(cpus)
		self._axeses = {}

	def subplot(self, *args):
		"""
		Create a subplot or return the axes previously created with the same
		args tuple. Either way, make it the current axes.

		The caller must use consistent args for each subplot, e.g. `(2, 2, 1)`
		vs. `(221,)`. Don't intermix them.

		This does not yet support keyword args.

		Use this to work around this matplotlib 2.2.2 deprecation when you want
		to retrieve a previously-created axes (subplot):

			MatplotlibDeprecationWarning: Adding an axes using the same
			arguments as a previous axes currently reuses the earlier instance.
			In a future version, a new instance will always be created and
			returned.  Meanwhile, this warning can be suppressed, and the
			future behavior ensured, by passing a unique label to each axes
			instance.

		Args:
			*args (Tuple[int]): subplot grid and index args

		Returns:
			(axes): the created or retrieved axes object; now the current one
		"""
		axes = self._axeses.get(args)
		if not axes:
			self._axeses[args] = axes = plt.subplot(*args)
		plt.sca(axes)
		return axes

	@staticmethod
	def set_ylim(axes, ymin, ymax):
		"""
		Set the axes y-limits, avoiding the matplotlib warning:

			UserWarning: Attempting to set identical bottom==top results

		Args:
			axes (axes): the axes to modify
			ymin (float): the minimum y value
			ymax (float): the maximum y value

		Returns:
			ylimits (Tuple[float]): the new y-axis limits as (`bottom`, `top`)
		"""
		if ymin == ymax:
			ymin -= 0.001
			ymax += 0.001
		return axes.set_ylim(ymin, ymax)

	@abc.abstractmethod
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		"""Inner method that each analysis class must override."""
		raise NotImplementedError("AnalysisPlot subclass must implement do_plot()")

	def plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		"""Public method to set up, make a plot, and cleanup."""
		fp.makedirs(plotOutDir)  # TODO(jerry): don't repeat this in 132 do_plot() methods

		with memory_debug.detect_leaks(), mp.rc_context():
			self.do_plot(inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata)

		self._axeses = {}

	@classmethod
	def main(cls, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile=None, metadata=None, cpus=0):
		"""Run an analysis plot for a Firetask."""
		instance = cls(cpus)
		instance.plot(inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata)
