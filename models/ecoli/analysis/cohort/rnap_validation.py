"""
Plots the RNA Polymerase abundances for all seeds and all generations as a
violin plot.

Notes
-----
There are hard-coded validation values:
	- one represents the average cell (defined as 44% along the cell cycle in
	  age - Neidhardt et al. Physiology of the Bacterial Cell. 1931). To be
	  comparable with the validation data, the molecule abundances are taken
	  from 44% along the cell cycle.
	- RNA polymerase abundance reported in Bremer & Dennis 2008, Table 3, for
	  the 40 minute and 60 minute doubling cells.
"""

from __future__ import absolute_import

import os
import cPickle
import numpy as np
from scipy import interpolate
from multiprocessing import Pool
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.utils import units, parallelization

from models.ecoli.analysis import cohortAnalysisPlot

# First generation (counting from zero) from which to gather doubling time
# values.  If fewer generations were run, this script quits early without
# plotting anything.
FIRST_GENERATION = 2

CELL_CYCLE_FRACTION = 0.44 # Average cell is 44% along its cell cycle length
RNAP_VALIDATION = {        # Bremer & Dennis 2008, Table 3
	'doubling_time': [20, 24, 30, 40, 60, 100],
    'rnap_abundance': [10.2e3, 10.0e3, 8.4e3, 5.7e3, 3.5e3, 1.8e3]}


FIGSIZE = (2.5, 5)
COUNTS_BOUNDS = [0, 7e3]


def mp_worker(sim_dir):
	sim_out_dir = os.path.join(sim_dir, 'simOut')
	rnap_count_avg_cell = None

	try:
		bulk_molecule_reader = TableReader(os.path.join(sim_out_dir, 'BulkMolecules'))
		index_rnap = bulk_molecule_reader.readAttribute('objectNames').index(rnap_id)
		rnap_count = bulk_molecule_reader.readColumn('counts', np.array([index_rnap]))

		unique_molecule_reader = TableReader(os.path.join(sim_out_dir, 'UniqueMoleculeCounts'))
		unique_molecule_ids = unique_molecule_reader.readAttribute('uniqueMoleculeIds')
		unique_molecule_counts = unique_molecule_reader.readColumn('uniqueMoleculeCounts')
		unique_molecule_reader.close()

		index_rnap = unique_molecule_ids.index('activeRnaPoly')
		rnap_active_count = unique_molecule_counts[:, index_rnap]

		index_average_cell = int(len(rnap_active_count) * CELL_CYCLE_FRACTION)
		rnap_count_avg_cell = rnap_count[index_average_cell] + rnap_active_count[index_average_cell]

	except Exception as e:
		print('Excluded from analysis due to broken files: {}'.format(sim_out_dir))

	return rnap_count_avg_cell


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(variantDir):
			raise Exception, 'variantDir does not currently exist as a directory'

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		analysis_paths = AnalysisPaths(variantDir, cohort_plot = True)
		n_gens = analysis_paths.n_generation

		# Check for sufficient generations
		if n_gens - 1 < FIRST_GENERATION:
			print 'Not enough generations to plot.'
			return

		sim_dirs = analysis_paths.get_cells(
			generation=range(FIRST_GENERATION, n_gens), seed = range(8))

		sim_data = cPickle.load(open(simDataFile, 'rb'))

		global rnap_id
		rnap_id = sim_data.moleculeIds.rnapFull

		p = Pool(parallelization.cpus())
		output = p.map(mp_worker, sim_dirs)
		p.close()
		p.join()

		# Filter output from broken files
		rnap_counts = [x for x in output if x]

		if not len(rnap_counts):
			print('Skipping plot due to no viable sims.')
			return

		# Plot
		doubling_time = sim_data.conditionToDoublingTime[sim_data.condition].asNumber(units.min)
		params = interpolate.splrep(
			RNAP_VALIDATION['doubling_time'],
			RNAP_VALIDATION['rnap_abundance'])
		rnap_abundance_fit = interpolate.splev(doubling_time, params)

		fig, ax = plt.subplots(1, 1, figsize=FIGSIZE)
		ax.violinplot(rnap_counts)
		ax.axhline(rnap_abundance_fit, color='tab:orange', lw=1)
		ax.set_ylim(*COUNTS_BOUNDS)
		ax.set_xlim([0.5, 1.5])
		ax.set_xticks([])
		y_ticks = ax.get_yticks()
		ax.set_yticklabels([])
		ax.spines['right'].set_visible(False)
		exportFigure(plt, plotOutDir, '{}__clean'.format(plotOutFileName), None)

		ax.set_title('n = {}'.format(len(rnap_counts)))
		ax.set_ylabel('Molecule abundance (counts)')
		ax.set_yticks(y_ticks)
		ax.spines['right'].set_visible(True)
		ax.spines['left'].set_visible(True)
		plt.subplots_adjust(left=0.4, bottom=0.2, right=0.6, top=0.8)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()