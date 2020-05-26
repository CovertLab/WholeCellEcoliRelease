"""
Common code for Multigen analysis plots.
"""

from __future__ import absolute_import, division, print_function

from models.ecoli.analysis import analysisPlot
from runscripts.manual import analysisMultigen


# noinspection PyAbstractClass
class MultigenAnalysisPlot(analysisPlot.AnalysisPlot):
	"""Abstract base class for Multigen analysis plots."""

	def cli(self):
		"""Run the command line for this AnalysisPlot subclass."""
		script = analysisMultigen.AnalysisMultigen(analysis_plotter=self)
		script.cli()
