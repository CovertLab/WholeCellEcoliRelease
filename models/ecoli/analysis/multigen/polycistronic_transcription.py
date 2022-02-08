"""
Plots the dynamics of the expression of a chosen polycistronic transcript.
"""

import pickle
import os

from matplotlib import pyplot as plt

from models.ecoli.analysis import multigenAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure, read_stacked_columns)
from wholecell.io.tablereader import TableReader

TU_ID = 'TU0-13388[c]'


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		if metadata['operons'] == 'off':
			print('Skipping analysis - no operons used in this sim')
			return

		cell_paths = self.ap.get_cells()

		simOutDir = os.path.join(cell_paths[0], "simOut")
		mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'mRNACounts'))
		mRNA_ids = mRNA_counts_reader.readAttribute('mRNA_ids')
		mRNA_cistron_ids = mRNA_counts_reader.readAttribute('mRNA_cistron_ids')
		monomer_counts_reader = TableReader(os.path.join(simOutDir, 'MonomerCounts'))
		all_monomer_ids = monomer_counts_reader.readAttribute('monomerIds')

		# Get indexes of cistrons that constitute the TU
		tu_index = mRNA_ids.index(TU_ID)
		cistron_ids = [
			sim_data.process.transcription.cistron_data['id'][i]
			for i in sim_data.process.transcription.rna_id_to_cistron_indexes(TU_ID)]
		cistron_indexes = [
			mRNA_cistron_ids.index(cistron_id) for cistron_id in cistron_ids]

		# Get indexes of monomers translated from the TU
		cistron_id_to_monomer_id = {
			monomer['cistron_id']: monomer['id']
			for monomer in sim_data.process.translation.monomer_data
			}
		monomer_ids = [
			cistron_id_to_monomer_id[cistron_id] for cistron_id in cistron_ids]
		monomer_indexes = [
			all_monomer_ids.index(monomer_id) for monomer_id in monomer_ids]

		# Load data
		time = read_stacked_columns(cell_paths, 'Main', 'time')
		tu_counts = read_stacked_columns(
			cell_paths, 'mRNACounts', 'mRNA_counts')[:, tu_index]
		cistron_counts = read_stacked_columns(
			cell_paths, 'mRNACounts', 'mRNA_cistron_counts')[:, cistron_indexes]
		monomer_counts = read_stacked_columns(
			cell_paths, 'MonomerCounts', 'monomerCounts')[:, monomer_indexes]

		plt.figure(figsize=(8.5, 11))

		# Plot counts of transcription unit
		ax = plt.subplot(3, 1, 1)
		ax.plot(time/60, tu_counts, label=TU_ID)
		ax.legend()
		ax.set_ylabel('Counts')

		# Plot counts of each cistron
		ax = plt.subplot(3, 1, 2)
		for i, cistron_id in enumerate(cistron_ids):
			ax.plot(time/60, cistron_counts[:, i], label=cistron_id)
		ax.legend()
		ax.set_ylabel('Counts')

		# Plot counts of each monomer
		ax = plt.subplot(3, 1, 3)
		for i, monomer_id in enumerate(monomer_ids):
			ax.plot(time/60, monomer_counts[:, i], label=monomer_id)
		ax.legend()
		ax.set_ylabel('Counts')
		ax.set_xlabel('Time (min)')

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
