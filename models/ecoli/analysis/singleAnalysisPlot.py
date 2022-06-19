"""
Common code for Single analysis plots.
"""

from __future__ import absolute_import, division, print_function

from models.ecoli.analysis import analysisPlot
from runscripts.manual import analysisSingle


# noinspection PyAbstractClass
class SingleAnalysisPlot(analysisPlot.AnalysisPlot):
	"""Abstract base class for Single analysis plots."""

	def cli(self):
		"""Run the command line for this AnalysisPlot subclass."""
		script = analysisSingle.AnalysisSingle(analysis_plotter=self)
		script.cli()
