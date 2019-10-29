"""
Comparison of metabolite concentrations that are not produced without CdsA.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 9/23/19
"""

from __future__ import absolute_import, division, print_function

import os

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import filepath
from wholecell.utils.sparkline import whitePadSparklineAxis


LIMITED_METABOLITES = [
	'CPD-12819[c]',
	'CPD-8260[c]',
]


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(variantDir):
			raise Exception, 'variantDir does not currently exist as a directory'

		filepath.makedirs(plotOutDir)

		ap = AnalysisPaths(variantDir, cohort_plot=True)

		limited_metabolites = []
		for sim_dir in ap.get_cells():
			sim_out_dir = os.path.join(sim_dir, 'simOut')

			# Listeners used
			kinetics_reader = TableReader(os.path.join(sim_out_dir, "EnzymeKinetics"))

			# Load data
			try:
				metabolite_indices = {m: i for i, m in enumerate(kinetics_reader.readAttribute('metaboliteNames'))}
				metabolite_counts = kinetics_reader.readColumn("metaboliteCountsFinal")[1:, :]
				counts_to_molar = kinetics_reader.readColumn('countsToMolar')[1:].reshape(-1, 1)
			except:
				print('Error reading data from {}'.format(sim_out_dir))
				continue

			# Calculate concentrations
			met_idx = np.array([metabolite_indices[m] for m in LIMITED_METABOLITES])
			metabolite_conc = counts_to_molar * metabolite_counts[:, met_idx]
			limited_metabolites += [metabolite_conc]

		limited_metabolites = np.vstack(limited_metabolites)

		# Values to calculate significance between different cohorts
		print('Metabolites: {}'.format(LIMITED_METABOLITES))
		print('Means: {}'.format(limited_metabolites.mean(axis=0)))
		print('Stds: {}'.format(limited_metabolites.std(axis=0)))
		print('N: {}'.format(limited_metabolites.shape[0]))

		plt.figure(figsize=(4, 4))
		xticks = [0, 1]

		# Plot data
		plt.violinplot(limited_metabolites, xticks, showmeans=True)

		# Format axes
		plt.ylim([0, 50])
		whitePadSparklineAxis(plt.gca())
		plt.xticks(xticks, LIMITED_METABOLITES)
		plt.ylabel('Concentration (uM)')

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
