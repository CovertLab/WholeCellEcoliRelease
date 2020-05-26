"""
Common code for Variant analysis plots.
"""

from __future__ import absolute_import, division, print_function

from models.ecoli.analysis import analysisPlot
from runscripts.manual import analysisVariant


# noinspection PyAbstractClass
class VariantAnalysisPlot(analysisPlot.AnalysisPlot):
	"""Abstract base class for Variant analysis plots."""

	def cli(self):
		"""Run the command line for this AnalysisPlot subclass."""
		script = analysisVariant.AnalysisVariant(analysis_plotter=self)
		script.cli()
