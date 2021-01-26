"""
Plot of outer membrane major components expression
Note: This is not an exhaustive list of all components of outer membrane.
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

outer_mem_proteins = {
	# phospholipids
	'CPD-12819[c]': 'phosphatidylathanolamine',
    'CPD-12824[c]' : 'cardiolipin',
    'CPD-8260[c]': 'phosphatidylglycerol',

    'CPD0-939[c]': 'LPS',

	# outer membrane proteins
    'EG10669-MONOMER[i]': 'ompA',
    # 'EG10671-MONOMER[e]': 'ompF',
    # 'EG10670-MONOMER[e]': 'ompC',
    # 'G6657-MONOMER[o]': 'ompG',
    # 'G7814-MONOMER[o]': 'ompL',
    # 'G6700-MONOMER[o]': 'ompN',
    'EG10673-MONOMER[o]': 'ompT',
    'EG12117-MONOMER[o]': 'ompX',

	# outer membrane porins
	# 'CPLX0-7530[e]': 'outer membrane porin E',
	'CPLX0-7533[e]': 'outer membrane porin C',
	'CPLX0-7534[e]': 'outer membrane porin F',
	# 'CPLX0-7655[e]': 'maltose outer membrane porin / phage lambda receptor protein',

	# Tol pal system
	# 'CPLX0-2201[s]': 'Tol-Pal system', # spans envelope
	# 'EG10684-MONOMER[e]': 'Pal',

	# other outer membrane
	# 'CPLX0-7944[e]': 'outer membrane phospholipase',
	# 'CPLX0-7952[s]': 'ferric coprogen outer membrane transport complex',
	# 'CPLX0-7964[e]': 'TolC outer membrane channel',
	# 'CPLX0-8009[e]': 'outer membrane protein; export and assembly of type 1 fimbriae',
	# 'CPLX0-1924[s]': 'vitamin B12 outer membrane transport complex',
	# 'CPLX0-1941[s]': 'ferric enterobactin outer membrane transport complex',
	# 'CPLX0-1942[s]': 'ferrichrome outer membrane transport complex',
	# 'CPLX0-1943[s]': 'ferric citrate outer membrane transport complex',
	'EG10544-MONOMER[o]': 'murein lipoprotein',

	# # phospholipid transport
    # 'EG11293-MONOMER[o]': 'lolB',
    # 'G6465-MONOMER[p]': 'lolA',

	# LPS biosynthesis
	'EG11424-MONOMER[i]': 'waaL',
	'EG10613-MONOMER[i]': 'msbA',

	# LPS transport
	'CPLX0-7704[i]': 'LPS ABC transporter (lipid A-core flippase)',
	'CPLX0-7992[i]': 'LPS transport system',
	'ABC-53-CPLX[i]': 'lptB2-F-G',
	'CPLX0-7628[e]': 'lptE-D',
	'YHBN-MONOMER[e]': 'lptA',
	'YHBG-MONOMER[i]': 'lptB',
	'G7664-MONOMER[i]': 'lptC',
	'EG11569-MONOMER[e]': 'lptD',
	'EG10855-MONOMER[e]': 'lptE',
	'G7888-MONOMER[i]': 'lptF',
	'G7889-MONOMER[i]': 'lptG',
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
											  (outer_mem_proteins.keys(),))
		counts = counts.astype(float).T

		row_max = counts.max(axis=1)
		normalized_counts = counts / row_max[:, np.newaxis]

		fig, ax = plt.subplots()
		image = ax.imshow(normalized_counts,
			cmap='RdBu',
			interpolation='nearest',
			aspect='auto',
			extent=[time[0],time[-1],len(outer_mem_proteins)-0.5,-0.5])
		ax.set_yticks(np.arange(0, len(outer_mem_proteins), 1))
		ax.set_yticklabels([outer_mem_proteins[mol_id] for
							mol_id in outer_mem_proteins.keys()], fontsize=8)
		plt.xlabel('time (s)')
		plt.title('Outer membrane major components expression')

		# add second axes to display final counts
		final_counts = counts[:, -1].astype(int).tolist()
		final_counts.reverse()

		ax2 = fig.add_subplot(111, sharex=ax, frameon=False)
		ax2.set_ylim([-0.5, len(outer_mem_proteins)-0.5])
		ax2.set_yticks(np.arange(0, len(outer_mem_proteins), 1))
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
