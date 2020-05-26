"""
Common code for Cohort analysis plots.
"""

from __future__ import absolute_import, division, print_function

from models.ecoli.analysis import analysisPlot
from runscripts.manual import analysisCohort


# noinspection PyAbstractClass
class CohortAnalysisPlot(analysisPlot.AnalysisPlot):
	"""Abstract base class for Cohort analysis plots."""

	def cli(self):
		"""Run the command line for this AnalysisPlot subclass."""
		script = analysisCohort.AnalysisCohort(analysis_plotter=self)
		script.cli()
