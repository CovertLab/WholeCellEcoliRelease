"""
Plots a violin/swarm plot of the excess protein monomer index that quantifies
how much excess monomers are produced for each protein complex.
"""

from functools import reduce
import pickle
import os

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		# Load from sim_data
		all_subunit_ids = sim_data.process.complexation.molecule_names
		stoich_matrix_monomers = sim_data.process.complexation.stoich_matrix_monomers()
		monomer_id_to_cistron_id = {
			monomer['id']: monomer['cistron_id']
			for monomer in sim_data.process.translation.monomer_data
			}
		cistron_id_to_rna_indexes = sim_data.process.transcription.cistron_id_to_rna_indexes

		# Listeners used
		monomer_counts_reader =	TableReader(os.path.join(simOutDir, 'MonomerCounts'))

		# Load data
		all_monomer_counts = monomer_counts_reader.readColumn('monomerCounts')
		monomer_ids = monomer_counts_reader.readAttribute('monomerIds')
		monomer_id_to_index = {
			monomer_id: i for (i, monomer_id) in enumerate(monomer_ids)}

		# Take averages across timepoints
		all_monomer_counts_mean = all_monomer_counts.mean(axis=0)

		excess_monomer_index = []
		cotranscribed_mask = []

		def is_cotranscribed(subunit_ids):
			"""
			Returns True if the cistrons that each encode for the list of
			protein monomer subunits are ever found on the same transcription
			unit.
			"""
			# Get IDs of cistrons that encode for each monomer
			cistron_ids = [
				monomer_id_to_cistron_id[monomer_id]
				for monomer_id in subunit_ids]

			# Get indexes of transcription units that contain each cistron
			rna_indexes = [
				cistron_id_to_rna_indexes(cistron_id)
				for cistron_id in cistron_ids]

			# Find any overlapping transcription unit indexes
			return len(reduce(np.intersect1d, rna_indexes)) > 0

		# Loop through each protein complex
		for complex_index in range(stoich_matrix_monomers.shape[1]):
			# Get stoichiometries and IDs of each monomer subunit
			subunit_mask = stoich_matrix_monomers[:, complex_index] < 0
			subunit_stoichs = -stoich_matrix_monomers[:, complex_index][subunit_mask]
			subunit_indexes = np.where(subunit_mask)[0]
			subunit_ids = [all_subunit_ids[i] for i in subunit_indexes]

			# Skip homogeneous protein complexes
			if len(subunit_indexes) == 1:
				continue

			# Get indexes of monomer subunits in the monomer counts table
			try:
				monomer_indexes = np.array([
					monomer_id_to_index[subunit_id] for subunit_id in subunit_ids
					])
			except KeyError:
				# Skip complexes with non-protein monomers or monomers that
				# were excluded from the model
				continue

			# Get counts of monomers and rescale with stoichiometries
			monomer_counts = all_monomer_counts_mean[monomer_indexes]
			rescaled_monomer_counts = monomer_counts/subunit_stoichs

			# Skip complexes with zero counts
			if np.all(rescaled_monomer_counts == 0):
				continue

			# Calculate excess monomer index and determine if monomers are
			# cotranscribed
			excess_monomer_index.append(
				(rescaled_monomer_counts.max() - rescaled_monomer_counts.min())
				/ rescaled_monomer_counts.sum()
				)
			cotranscribed_mask.append(
				is_cotranscribed(subunit_ids))

		excess_monomer_index = np.array(excess_monomer_index)
		cotranscribed_mask = np.array(cotranscribed_mask)

		# Overlay swarm and violin plots
		plt.figure(figsize=(3, 6))
		x_jitter = np.random.normal(1, 0.01, len(excess_monomer_index))
		plt.scatter(
			x_jitter[~cotranscribed_mask],
			excess_monomer_index[~cotranscribed_mask],
			s=5, label='not cotranscribed')
		plt.scatter(
			x_jitter[cotranscribed_mask],
			excess_monomer_index[cotranscribed_mask],
			s=5, label='cotranscribed')

		plt.violinplot(excess_monomer_index[~cotranscribed_mask])
		if np.any(cotranscribed_mask):
			plt.violinplot(excess_monomer_index[cotranscribed_mask])

		plt.legend()
		plt.xticks([])
		plt.ylabel('Excess monomer index')
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
