"""
Common code for causality network analysis scripts.
"""

from __future__ import absolute_import
from __future__ import division

from models.ecoli.analysis import analysisPlot
from runscripts.manual import buildCausalityNetwork


class CausalityNetworkAnalysis(analysisPlot.AnalysisPlot):
	"""Abstract base class for causality network analysis scripts."""

	def cli(self):
		"""Run the command line for this AnalysisPlot subclass."""
		script = buildCausalityNetwork.BuildCausalityNetwork(analysis_plotter=self)
		script.cli()
