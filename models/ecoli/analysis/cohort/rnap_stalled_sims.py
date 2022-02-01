"""
Cohort analysis identifying seeds and generations in which
stalled transcript elongation occurred
"""

import os
import numpy as np
from models.ecoli.analysis import cohortAnalysisPlot
from wholecell.io.tablereader import TableReader


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		cell_paths = self.ap.get_cells()

		seed_gen_stalled = []
		for sim_dir in cell_paths:
			simOutDir = os.path.join(sim_dir, 'simOut')

			# Listeners used
			main_reader = TableReader(os.path.join(simOutDir, 'Main'))
			rnap_data_reader = TableReader(os.path.join(simOutDir, 'RnapData'))

			# Load data
			time = main_reader.readColumn('time')
			stalled_rnap_counts = rnap_data_reader.readColumn('didStall')

			# add to result if stalled elongation occurred
			if np.sum(stalled_rnap_counts) > 0:
				paths = os.path.normpath(simOutDir).split(os.sep)
				seed = paths[-4]
				gen = paths[-3].split('_')[1]
				seed_gen_stalled.append([seed, gen, simOutDir])

		seed_gen_stalled.sort()
		np.savetxt(fname=f"{plotOutDir}/rnap_stalled_sims.tsv",
				   X=seed_gen_stalled, delimiter='\t', fmt='%s',
				   comments='', header='seed\tgeneration\tsimOutDir')


if __name__ == '__main__':
	Plot().cli()
