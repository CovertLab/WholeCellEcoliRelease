"""
Common code for causality network analysis scripts.
"""

from __future__ import absolute_import, division, print_function

from models.ecoli.analysis import analysisPlot
from runscripts.manual import buildCausalityNetwork


# noinspection PyAbstractClass
class CausalityNetworkAnalysis(analysisPlot.AnalysisPlot):
	"""Abstract base class for causality network analysis scripts."""

	def cli(self):
		"""Run the command line for this AnalysisPlot subclass."""
		script = buildCausalityNetwork.BuildCausalityNetwork(analysis_plotter=self)
		script.cli()
