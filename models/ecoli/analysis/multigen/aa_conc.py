"""
Time series and histogram of amino acid concentrations

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/17/19
"""

from __future__ import absolute_import
from __future__ import division

import cPickle
import os

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import multigenAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.analysis.analysis_tools import read_bulk_molecule_counts
from wholecell.io.tablereader import TableReader
from wholecell.utils import filepath
from wholecell.utils import units


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(seedOutDir):
			raise Exception, 'seedOutDir does not currently exist as a directory'

		filepath.makedirs(plotOutDir)

		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)

		aa_ids = sim_data.moleculeGroups.aaIDs
		targets = np.array([sim_data.process.metabolism.concDict[key].asNumber(units.mmol / units.L) for key in aa_ids])

		ap = AnalysisPaths(seedOutDir, multi_gen_plot=True)

		aa_conc = None
		time = None
		for sim_dir in ap.get_cells():
			simOutDir = os.path.join(sim_dir, 'simOut')

			# Load data
			## Readers
			main_reader = TableReader(os.path.join(simOutDir, 'Main'))
			kinetics_reader = TableReader(os.path.join(simOutDir, 'EnzymeKinetics'))

			## Read data
			(aa_counts,) = read_bulk_molecule_counts(simOutDir, aa_ids)
			counts_to_molar = kinetics_reader.readColumn('countsToMolar').reshape(-1, 1)
			cell_time = main_reader.readColumn('time')

			## Append for all generations
			if aa_conc is None:
				aa_conc = aa_counts[1:, :]*counts_to_molar[1:]
				time = cell_time[1:] / 60
			else:
				aa_conc = np.vstack((aa_conc, aa_counts[1:, :]*counts_to_molar[1:]))
				time = np.hstack((time, cell_time[1:] / 60))

		# Plot time series
		plt.figure(figsize=(8.5, 11))

		for idx in xrange(21):
			plt.subplot(6, 4, idx + 1)
			plt.plot(time, aa_conc[:, idx])
			plt.ylim(0)
			plt.axhline(targets[idx], color='r', linestyle='--')
			plt.title(aa_ids[idx], fontsize=8)
			plt.tick_params(labelsize=8)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		# Plot histogram
		plt.figure(figsize=(8.5, 11))

		n_row = 6
		n_col = 4
		for idx in xrange(len(aa_ids)):
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
