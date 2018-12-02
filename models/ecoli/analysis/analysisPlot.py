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

from __future__ import absolute_import
from __future__ import division

import abc
import matplotlib as mp
from wholecell.utils import memory_debug, parallelization


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
		self.cpus = parallelization.cpus(cpus)

	@abc.abstractmethod
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		"""Inner method that each analysis class must override."""
		raise NotImplementedError("AnalysisPlot subclass must implement do_plot()")

	def plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		"""Public method to set up, make a plot, and cleanup."""
		with memory_debug.detect_leaks(), mp.rc_context():
			self.do_plot(inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata)

	@classmethod
	def main(cls, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile=None, metadata=None, cpus=0):
		"""Run an analysis plot for a Firetask."""
		instance = cls(cpus)
		instance.plot(inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata)
