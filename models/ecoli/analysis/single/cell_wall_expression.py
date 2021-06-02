"""
Plot of cell wall major components expression
Note: This is not an exhaustive list of all components of cell wall.
"""

import pickle
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

CELL_WALL_PROTEINS = {
    # murein
    'CPD-12261': 'murein',
    'CPD-12231':'peptidoglycan dimer',
    'C6':'lipid II',
    'C5': 'lipid I',
    'CPD-9646': 'undecaprenyl phosphate',
    'UNDECAPRENYL-DIPHOSPHATE': 'undecaprenyl diphosphate',

    # catalyzes lipid II transfer to outer membrane
    'G6561-MONOMER': 'lipid II flippase (MurJ)',
    'EG10344-MONOMER': 'essential cell division protein (FtsW)',

    # catalyzes lipid II formation
    'NACGLCTRANS-MONOMER': 'MurG',

    # catalyzes lipid I formation
    'PHOSNACMURPENTATRANS-MONOMER': 'MraY',

    # PBP
    'CPLX0-7717': 'PBP1A', # transglycosylase-transpeptidase ~100
    'CPLX0-3951': 'PBP1B', # transglycosylase-transpeptidase ~100
    'G7322-MONOMER': 'PBP1C', # transglycosylase
    'EG10606-MONOMER': 'PBP2', # transpeptidase ~20
    'EG10341-MONOMER': 'PBP3', # transglycosylase-transpeptidase ~50
    'EG10202-MONOMER': 'PBP4', # DD-endopeptidase, DD-carboxypeptidase ~110
    'EG10201-MONOMER': 'PBP5', # DD-caroxypeptidase ~1,800
    'EG10203-MONOMER': 'PBP6', # DD-carbocypeptidase ~600
    'EG12015-MONOMER': 'PBP7', # DD-endopeptidase

}

class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	_suppress_numpy_warnings = True

	def do_plot(self, simOutDir, plotOutDir, plotOutFileName,
				simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		# Listeners used
		main_reader = TableReader(os.path.join(simOutDir, 'Main'))

		# Load data
		initial_time = main_reader.readAttribute('initialTime')
		time = main_reader.readColumn('time') - initial_time

		existing_molecule_ids = [
			mol_id + sim_data.getter.get_compartment_tag(mol_id)
			for mol_id in CELL_WALL_PROTEINS.keys()
			if sim_data.getter.is_valid_molecule(mol_id)
		]
		(counts,) = read_bulk_molecule_counts(
			simOutDir, (existing_molecule_ids,))
		counts = counts.astype(float).T

		row_max = counts.max(axis=1)
		normalized_counts = counts / row_max[:, np.newaxis]

		fig, ax = plt.subplots()
		image = ax.imshow(normalized_counts,
			cmap='RdBu',
			interpolation='nearest',
			aspect='auto',
			extent=[time[0],time[-1],len(existing_molecule_ids)-0.5,-0.5])
		ax.set_yticks(np.arange(0, len(existing_molecule_ids), 1))
		ax.set_yticklabels([CELL_WALL_PROTEINS[mol_id[:-3]] for
							mol_id in existing_molecule_ids], fontsize=8)
		plt.xlabel('time (s)')
		plt.title('Cell wall major components expression')

		# add second axes to display final counts
		final_counts = counts[:, -1].astype(int).tolist()
		final_counts.reverse()

		ax2 = fig.add_subplot(111, sharex=ax, frameon=False)
		ax2.set_ylim([-0.5, len(existing_molecule_ids)-0.5])
		ax2.set_yticks(np.arange(0, len(existing_molecule_ids), 1))
		ax2.set_yticklabels(final_counts, fontsize=8)
		ax2.tick_params(length=0)
		ax2.yaxis.tick_right()
		ax2.yaxis.set_label_position("right")
		ax2.get_xaxis().set_visible(False)
		ax2.set_ylabel('final counts')

		ax3 = fig.add_axes([0.93, 0.12, 0.02, 0.28])
		plt.colorbar(image, cax=ax3)

		plt.subplots_adjust(left=0.4, right=0.75)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
