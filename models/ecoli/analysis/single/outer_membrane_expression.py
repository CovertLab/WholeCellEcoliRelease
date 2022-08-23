"""
Plot of outer membrane major components expression
Note: This is not an exhaustive list of all components of outer membrane.
"""

from __future__ import absolute_import, division, print_function

import os
import pickle

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.analysis.analysis_tools import read_bulk_molecule_counts
from wholecell.io.tablereader import TableReader
from wholecell.analysis.plotting_tools import CMAP_COLORS_255

CMAP_COLORS = [[shade/255. for shade in color] for color in CMAP_COLORS_255]
CMAP_OVER = [0, 1, 0.75]

OUTER_MEM_PROTEINS = {
	# phospholipids
	'CPD-12819': 'phosphatidylathanolamine',
    'CPD-12824' : 'cardiolipin',
    'CPD-8260': 'phosphatidylglycerol',
    'CPD0-939': 'LPS',

	# outer membrane proteins
    'EG10669-MONOMER': 'ompA',
    'EG10671-MONOMER': 'ompF',
    'EG10670-MONOMER': 'ompC',
    'G6657-MONOMER': 'ompG',
    'G7814-MONOMER': 'ompL',
    'G6700-MONOMER': 'ompN',
    'EG10673-MONOMER': 'ompT',
    'EG12117-MONOMER': 'ompX',

	# outer membrane porins
	'CPLX0-7530': 'outer membrane porin E',
	'CPLX0-7533': 'outer membrane porin C',
	'CPLX0-7534': 'outer membrane porin F',
	'CPLX0-7655': 'maltose outer membrane porin / phage lambda receptor protein',

	# Tol pal system
	'CPLX0-2201': 'Tol-Pal system', # spans envelope
	'EG10684-MONOMER': 'Pal',

	# other outer membrane
	'CPLX0-7944': 'outer membrane phospholipase',
	'CPLX0-7952': 'ferric coprogen outer membrane transport complex',
	'CPLX0-7964': 'TolC outer membrane channel',
	'CPLX0-8009': 'outer membrane protein; export and assembly of type 1 fimbriae',
	'CPLX0-1924': 'vitamin B12 outer membrane transport complex',
	'CPLX0-1941': 'ferric enterobactin outer membrane transport complex',
	'CPLX0-1942': 'ferrichrome outer membrane transport complex',
	'CPLX0-1943': 'ferric citrate outer membrane transport complex',
	'EG10544-MONOMER': 'murein lipoprotein',

	# # phospholipid transport
    'EG11293-MONOMER': 'lolB',
    'G6465-MONOMER': 'lolA',

	# LPS biosynthesis
	'EG11424-MONOMER': 'waaL',
	'EG10613-MONOMER': 'msbA',

	# LPS transport
	'CPLX0-7704': 'LPS ABC transporter (lipid A-core flippase)',
	'CPLX0-7992': 'LPS transport system',
	'ABC-53-CPLX': 'lptB2-F-G',
	'CPLX0-7628': 'lptE-D',
	'YHBN-MONOMER': 'lptA',
	'YHBG-MONOMER': 'lptB',
	'G7664-MONOMER': 'lptC',
	'EG11569-MONOMER': 'lptD',
	'EG10855-MONOMER': 'lptE',
	'G7888-MONOMER': 'lptF',
	'G7889-MONOMER': 'lptG',
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
			for mol_id in OUTER_MEM_PROTEINS.keys()
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
		ax.set_yticklabels([OUTER_MEM_PROTEINS[mol_id[:-3]] for
							mol_id in existing_molecule_ids], fontsize=8)
		plt.xlabel('time (s)')
		plt.title('Outer membrane major components expression')

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
