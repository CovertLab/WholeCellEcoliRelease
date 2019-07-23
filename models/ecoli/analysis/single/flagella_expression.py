"""
Plot of flagella protein expression

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 8/2/18
"""

from __future__ import absolute_import
from __future__ import division

import cPickle
import os

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.analysis.analysis_tools import read_bulk_molecule_counts
from wholecell.io.tablereader import TableReader
from wholecell.utils import filepath
from wholecell.analysis.plotting_tools import CMAP_COLORS_255

CMAP_COLORS = [[shade/255. for shade in color] for color in CMAP_COLORS_255]
CMAP_OVER = [0, 1, 0.75]

ordered_flagella_protein_ids = [
	'G7028-MONOMER[i]',
	'G378-MONOMER[s]',
	'G377-MONOMER[i]',
	'G370-MONOMER[i]',
	'EG11977-MONOMER[s]',
	'EG11976-MONOMER[s]',
	'EG11975-MONOMER[s]',
	'EG11656-MONOMER[s]',
	'EG11224-MONOMER[s]',
	'CPLX0-7451[s]',
	'FLIN-FLAGELLAR-C-RING-SWITCH[s]',
	'FLIM-FLAGELLAR-C-RING-SWITCH[s]',
	'FLIG-FLAGELLAR-SWITCH-PROTEIN[s]',
	'CPLX0-7450[s]',
	'FLGH-FLAGELLAR-L-RING[s]',
	'MOTA-FLAGELLAR-MOTOR-STATOR-PROTEIN[i]',
	'MOTB-FLAGELLAR-MOTOR-STATOR-PROTEIN[i]',
	'FLGB-FLAGELLAR-MOTOR-ROD-PROTEIN[s]',
	'FLGC-FLAGELLAR-MOTOR-ROD-PROTEIN[s]',
	'FLGF-FLAGELLAR-MOTOR-ROD-PROTEIN[s]',
	'FLGG-FLAGELLAR-MOTOR-ROD-PROTEIN[s]',
	'FLGI-FLAGELLAR-P-RING[s]',
	'FLIF-FLAGELLAR-MS-RING[s]',
	'EG11346-MONOMER[s]',
	'EG10322-MONOMER[s]',
	'FLAGELLAR-MOTOR-COMPLEX[s]',
	'G361-MONOMER[s]',
	'EG11967-MONOMER[s]',
	'EG11545-MONOMER[s]',
	'EG10321-MONOMER[s]',
	'EG10841-MONOMER[s]',
	'CPLX0-7452[s]',
]

flagella_proteins = {
	'G7028-MONOMER[i]': 'FlhB',
	'G378-MONOMER[s]': 'FliJ',
	'G377-MONOMER[i]': 'FliI',
	'G370-MONOMER[i]': 'FlhA',
	'EG11977-MONOMER[s]': 'FliR',
	'EG11976-MONOMER[s]': 'FliQ',
	'EG11975-MONOMER[s]': 'FliP',
	'EG11656-MONOMER[s]': 'FliH',
	'EG11224-MONOMER[s]': 'FliO',
	'CPLX0-7451[s]': 'export apparatus',
	'FLIN-FLAGELLAR-C-RING-SWITCH[s]': 'FliN',
	'FLIM-FLAGELLAR-C-RING-SWITCH[s]': 'FliM',
	'FLIG-FLAGELLAR-SWITCH-PROTEIN[s]': 'FliG',
	'CPLX0-7450[s]': 'motor switch complex',
	'FLGH-FLAGELLAR-L-RING[s]': 'FlgH',
	'MOTA-FLAGELLAR-MOTOR-STATOR-PROTEIN[i]': 'MotA',
	'MOTB-FLAGELLAR-MOTOR-STATOR-PROTEIN[i]': 'MotB',
	'FLGB-FLAGELLAR-MOTOR-ROD-PROTEIN[s]': 'FlgB',
	'FLGC-FLAGELLAR-MOTOR-ROD-PROTEIN[s]': 'FlgC',
	'FLGF-FLAGELLAR-MOTOR-ROD-PROTEIN[s]': 'FlgF',
	'FLGG-FLAGELLAR-MOTOR-ROD-PROTEIN[s]': 'FlgG',
	'FLGI-FLAGELLAR-P-RING[s]': 'FlgI',
	'FLIF-FLAGELLAR-MS-RING[s]': 'FliF',
	'EG11346-MONOMER[s]': 'FliE',
	'EG10322-MONOMER[s]': 'FliL',
	'FLAGELLAR-MOTOR-COMPLEX[s]': 'motor complex',
	'G361-MONOMER[s]': 'FlgE',
	'EG11967-MONOMER[s]': 'FlgK',
	'EG11545-MONOMER[s]': 'FlgL',
	'EG10321-MONOMER[s]': 'FliC',
	'EG10841-MONOMER[s]': 'FliD',
	'CPLX0-7452[s]': 'Flagellum',
}

class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, 'simOutDir does not currently exist as a directory'

		filepath.makedirs(plotOutDir)

		# Listeners used
		main_reader = TableReader(os.path.join(simOutDir, 'Main'))

		# Load data
		initial_time = main_reader.readAttribute('initialTime')
		time = main_reader.readColumn('time') - initial_time
		(counts,) = read_bulk_molecule_counts(simOutDir, (ordered_flagella_protein_ids,))
		counts = counts.astype(float).T

		row_max = counts.max(axis=1)
		normalized_counts = counts / row_max[:, np.newaxis]

		fig, ax = plt.subplots()
		image = ax.imshow(normalized_counts,
			cmap='hot',
			interpolation='nearest',
			aspect='auto',
			extent=[time[0],time[-1],len(ordered_flagella_protein_ids)-0.5,-0.5])
		ax.set_yticks(np.arange(0, len(flagella_proteins), 1))
		ax.set_yticklabels([flagella_proteins[mol_id] for mol_id in ordered_flagella_protein_ids], fontsize=8)
		plt.xlabel('time (s)')

		# add second axes to display final counts
		final_counts = counts[:, -1].astype(int).tolist()
		final_counts.reverse()

		ax2 = fig.add_subplot(111, sharex=ax, frameon=False)
		ax2.set_ylim([-0.5, len(ordered_flagella_protein_ids)-0.5])
		ax2.set_yticks(np.arange(0, len(flagella_proteins), 1))
		ax2.set_yticklabels(final_counts, fontsize=8)
		ax2.tick_params(length=0)
		ax2.yaxis.tick_right()
		ax2.yaxis.set_label_position("right")
		ax2.get_xaxis().set_visible(False)
		ax2.set_ylabel('final counts')

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
