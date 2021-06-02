"""
Plot of flagella protein expression
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

FLAGELLA_PROTEINS = {
	'G7028-MONOMER': 'FlhB',
	'G378-MONOMER': 'FliJ',
	'G377-MONOMER': 'FliI',
	'G370-MONOMER': 'FlhA',
	'EG11977-MONOMER': 'FliR',
	'EG11976-MONOMER': 'FliQ',
	'EG11975-MONOMER': 'FliP',
	'EG11656-MONOMER': 'FliH',
	'EG11224-MONOMER': 'FliO',
	'CPLX0-7451': 'export apparatus',
	'FLIN-FLAGELLAR-C-RING-SWITCH': 'FliN',
	'FLIM-FLAGELLAR-C-RING-SWITCH': 'FliM',
	'FLIG-FLAGELLAR-SWITCH-PROTEIN': 'FliG',
	'CPLX0-7450': 'motor switch complex',
	'FLGH-FLAGELLAR-L-RING': 'FlgH',
	'MOTA-FLAGELLAR-MOTOR-STATOR-PROTEIN': 'MotA',
	'MOTB-FLAGELLAR-MOTOR-STATOR-PROTEIN': 'MotB',
	'FLGB-FLAGELLAR-MOTOR-ROD-PROTEIN': 'FlgB',
	'FLGC-FLAGELLAR-MOTOR-ROD-PROTEIN': 'FlgC',
	'FLGF-FLAGELLAR-MOTOR-ROD-PROTEIN': 'FlgF',
	'FLGG-FLAGELLAR-MOTOR-ROD-PROTEIN': 'FlgG',
	'FLGI-FLAGELLAR-P-RING': 'FlgI',
	'FLIF-FLAGELLAR-MS-RING': 'FliF',
	'EG11346-MONOMER': 'FliE',
	'EG10322-MONOMER': 'FliL',
	'FLAGELLAR-MOTOR-COMPLEX': 'motor complex',
	'G361-MONOMER': 'FlgE',
	'EG11967-MONOMER': 'FlgK',
	'EG11545-MONOMER': 'FlgL',
	'EG10321-MONOMER': 'FliC',
	'EG10841-MONOMER': 'FliD',
	'CPLX0-7452': 'Flagellum',
}

class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	_suppress_numpy_warnings = True

	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		# Listeners used
		main_reader = TableReader(os.path.join(simOutDir, 'Main'))

		# Load data
		initial_time = main_reader.readAttribute('initialTime')
		time = main_reader.readColumn('time') - initial_time

		existing_molecule_ids = [
			mol_id + sim_data.getter.get_compartment_tag(mol_id)
			for mol_id in FLAGELLA_PROTEINS.keys()
			if sim_data.getter.is_valid_molecule(mol_id)
		]
		(counts,) = read_bulk_molecule_counts(simOutDir, (existing_molecule_ids,))
		counts = counts.astype(float).T

		row_max = counts.max(axis=1)
		normalized_counts = counts / row_max[:, np.newaxis]

		fig, ax = plt.subplots()
		image = ax.imshow(normalized_counts,
			cmap='hot',
			interpolation='nearest',
			aspect='auto',
			extent=[time[0],time[-1],len(existing_molecule_ids)-0.5,-0.5])
		ax.set_yticks(np.arange(0, len(existing_molecule_ids), 1))
		ax.set_yticklabels([FLAGELLA_PROTEINS[mol_id[:-3]] for mol_id in existing_molecule_ids], fontsize=8)
		plt.xlabel('time (s)')

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

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
