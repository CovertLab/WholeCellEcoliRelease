"""
Common code for Parca analysis plots.
"""

from __future__ import absolute_import, division, print_function

from models.ecoli.analysis import analysisPlot
from runscripts.manual import analysisParca


# noinspection PyAbstractClass
class ParcaAnalysisPlot(analysisPlot.AnalysisPlot):
	"""Abstract base class for Parca analysis plots."""

	def cli(self):
		"""Run the command line for this AnalysisPlot subclass."""
		script = analysisParca.AnalysisParca(analysis_plotter=self)
		script.cli()
