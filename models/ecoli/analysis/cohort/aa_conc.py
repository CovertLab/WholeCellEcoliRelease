"""
Comparison of average amino acid concentrations to expected concentrations
"""

from __future__ import absolute_import, division, print_function

import os

from matplotlib import pyplot as plt
import numpy as np
from six.moves import cPickle, range

from models.ecoli.analysis import cohortAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.analysis.analysis_tools import read_bulk_molecule_counts
from wholecell.io.tablereader import TableReader
from wholecell.utils import units


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)

		aa_ids = sim_data.molecule_groups.amino_acids
		targets = np.array([sim_data.process.metabolism.conc_dict[key].asNumber(units.mmol / units.L) for key in aa_ids])


		aa_conc = None
		for sim_dir in self.ap.get_cells():
			simOutDir = os.path.join(sim_dir, 'simOut')

			# Load data
			## Readers
			kinetics_reader = TableReader(os.path.join(simOutDir, 'EnzymeKinetics'))

			## Read data
			(aa_counts,) = read_bulk_molecule_counts(simOutDir, aa_ids)
			counts_to_molar = kinetics_reader.readColumn('countsToMolar').reshape(-1, 1)

			## Append for all generations
			if aa_conc is None:
				aa_conc = aa_counts[1:, :]*counts_to_molar[1:]
			else:
				aa_conc = np.vstack((aa_conc, aa_counts[1:, :]*counts_to_molar[1:]))

		ave_conc = aa_conc.mean(axis=0)
		max_conc = np.max((ave_conc, targets))

		# Plot correlation to data
		plt.figure()

		plt.loglog(targets, ave_conc, 'o')
		plt.loglog([0, max_conc], [0, max_conc], 'k--')

		plt.xlabel('Expected Amino Acid Conc (mM)')
		plt.ylabel('Average Simulation Conc (mM)')

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		# Plot histogram
		plt.figure(figsize=(8.5, 11))

		n_row = 6
		n_col = 4
		for idx in range(len(aa_ids)):
			ave_conc = aa_conc[:, idx].mean()
			plt.subplot(n_row, n_col, idx + 1)
			plt.hist(aa_conc[:, idx])
			plt.axvline(targets[idx], color='r', linestyle='--')
			plt.axvline(ave_conc, color='k', linestyle='--')
			plt.xlabel('Conc (mM)', fontsize=6)
			plt.title('{}\nAve: {:.3f} mM'.format(aa_ids[idx], ave_conc), fontsize=8)
			plt.tick_params(labelsize=8)

		plt.subplot(n_row, n_col, len(aa_ids) + 1)
		plt.axis('off')
		plt.plot(0, 0, 'r--')
		plt.plot(0, 0, 'k--')
		plt.legend(['Expected Conc', 'Average Simulation Conc'], fontsize=8, loc=10)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + '_hist', metadata)

		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
