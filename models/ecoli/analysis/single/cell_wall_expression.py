"""
Plot of cell wall major components expression
Note: This is not an exhaustive list of all components of cell wall.
"""

from __future__ import absolute_import, division, print_function

import os

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.analysis.analysis_tools import read_bulk_molecule_counts
from wholecell.io.tablereader import TableReader
from wholecell.analysis.plotting_tools import CMAP_COLORS_255


CMAP_COLORS = [[shade/255. for shade in color] for color in CMAP_COLORS_255]
CMAP_OVER = [0, 1, 0.75]

cell_wall_proteins = {
    # murein
    'CPD-12261[p]': 'murein',
    'CPD-12231[p]':'peptidoglycan dimer',
    'C6[p]':'lipid II',
    'C5[c]': 'lipid I',
    'CPD-9646[c]': 'undecaprenyl phosphate',
    'UNDECAPRENYL-DIPHOSPHATE[c]': 'undecaprenyl diphosphate',

    # catalyzes lipid II transfer to outer membrane
    'G6561-MONOMER[i]': 'lipid II flippase (MurJ)',
    'EG10344-MONOMER[i]': 'essential cell division protein (FtsW)',

    # catalyzes lipid II formation
    'NACGLCTRANS-MONOMER[c]': 'MurG',

    # catalyzes lipid I formation
    'PHOSNACMURPENTATRANS-MONOMER[i]': 'MraY',

    # PBP
    'CPLX0-7717[i]': 'PBP1A', # transglycosylase-transpeptidase ~100
    'CPLX0-3951[i]': 'PBP1B', # transglycosylase-transpeptidase ~100
    'G7322-MONOMER[i]': 'PBP1C', # transglycosylase
    'EG10606-MONOMER[i]': 'PBP2', # transpeptidase ~20
    'EG10341-MONOMER[i]': 'PBP3', # transglycosylase-transpeptidase ~50
    'EG10202-MONOMER[p]': 'PBP4', # DD-endopeptidase, DD-carboxypeptidase ~110
    'EG10201-MONOMER[i]': 'PBP5', # DD-caroxypeptidase ~1,800
    'EG10203-MONOMER[i]': 'PBP6', # DD-carbocypeptidase ~600
    'EG12015-MONOMER[p]': 'PBP7', # DD-endopeptidase

}

class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName,
				simDataFile, validationDataFile, metadata):
		# Listeners used
		main_reader = TableReader(os.path.join(simOutDir, 'Main'))

		# Load data
		initial_time = main_reader.readAttribute('initialTime')
		time = main_reader.readColumn('time') - initial_time
		(counts,) = read_bulk_molecule_counts(simOutDir,
											  (cell_wall_proteins.keys(),))
		counts = counts.astype(float).T

		row_max = counts.max(axis=1)
		normalized_counts = counts / row_max[:, np.newaxis]

		fig, ax = plt.subplots()
		image = ax.imshow(normalized_counts,
			cmap='RdBu',
			interpolation='nearest',
			aspect='auto',
			extent=[time[0],time[-1],len(cell_wall_proteins)-0.5,-0.5])
		ax.set_yticks(np.arange(0, len(cell_wall_proteins), 1))
		ax.set_yticklabels([cell_wall_proteins[mol_id] for
							mol_id in cell_wall_proteins.keys()], fontsize=8)
		plt.xlabel('time (s)')
		plt.title('Cell wall major components expression')

		# add second axes to display final counts
		final_counts = counts[:, -1].astype(int).tolist()
		final_counts.reverse()

		ax2 = fig.add_subplot(111, sharex=ax, frameon=False)
		ax2.set_ylim([-0.5, len(cell_wall_proteins.keys())-0.5])
		ax2.set_yticks(np.arange(0, len(cell_wall_proteins), 1))
		ax2.set_yticklabels(final_counts, fontsize=8)
		ax2.tick_params(length=0)
		ax2.yaxis.tick_right()
		ax2.yaxis.set_label_position("right")
		ax2.get_xaxis().set_visible(False)
		ax2.set_ylabel('final counts')

		ax3 = fig.add_axes([0.93, 0.12, 0.02, 0.28])
		plt.colorbar(image, cax=ax3)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
