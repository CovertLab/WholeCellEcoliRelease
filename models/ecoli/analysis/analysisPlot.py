"""
Common code for analysis plots. The abstract base class AnalysisPlot defines a
plot() method for scripts to call.

TODO: Set the fallback for the matplotlib back end so we needn't hack its
installation. Set the CWD so matplotlib finds the matplotlibrc file even when
FireWorks doesn't launch it in the right directory.

TODO: Setup/reset matplotlib before each script and cleanup afterwards.

TODO: Enable future warnings, esp. for matplotlib.

TODO: Move the run_plot() args to instance variables?

TODO: Other shared code to simplify the subclasses, e.g. make plotOutDir,
check that `os.path.isdir(simOutDir)`, instantiate an AnalysisPaths (except for
SingleAnalysisPlot subclasses), etc.
"""

from __future__ import absolute_import
from __future__ import division

import abc
from wholecell.utils import memory_debug


class AnalysisPlot(object):
	"""Abstract Base Class for analysis plots.

	Each analysis class must override do_plot().

	Call main() to run an analysis plot for a Firetask.

	Use the environment variable 'DEBUG_GC' to enable memory leak debugging.
	"""
	__metaclass__ = abc.ABCMeta

	@abc.abstractmethod
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		"""Inner method that each analysis class must override."""
		raise NotImplementedError("AnalysisPlot subclass must implement do_plot()")

	def plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		"""Public method to set up, make a plot, and cleanup."""
		with memory_debug.detect_leaks():
			# TODO: Setup.

			self.do_plot(inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata)

			# TODO: Cleanup.

	@classmethod
	def main(cls, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile=None, metadata=None):
		"""Run an analysis plot for a Firetask."""
		instance = cls()
		instance.plot(inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata)
