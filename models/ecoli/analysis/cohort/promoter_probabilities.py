"""
Compare expected promoter binding probabilities with observed probabilities.
Useful for checking if actual promoter binding differs from the expected and
if this is due to limitations from the counts of active TFs available to bind.
"""

import os
import pickle

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from wholecell.io.tablereader import TableReader


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	_suppress_numpy_warnings = True

	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		tf_to_gene_id = sim_data.process.transcription_regulation.tf_to_gene_id

		ap = AnalysisPaths(variantDir, cohort_plot=True)

		expected = []
		bound = []
		active_tfs = []
		promoters = []
		for sim_dir in ap.get_cells():
			simOutDir = os.path.join(sim_dir, 'simOut')

			# Listeners used
			rna_synth_reader = TableReader(os.path.join(simOutDir, 'RnaSynthProb'))

			# Load data
			expected.append(rna_synth_reader.readColumn('pPromoterBound'))
			bound.append(rna_synth_reader.readColumn('nActualBound'))
			promoters.append(rna_synth_reader.readColumn('n_available_promoters'))

			tf_ids = rna_synth_reader.readAttribute('tf_ids')
			tf_ids_with_location = [f'{tf}[c]' for tf in tf_ids]
			active_counts, = read_bulk_molecule_counts(simOutDir, (tf_ids_with_location,))
			active_tfs.append(active_counts)

		# Turn lists into array for easier processing
		expected = np.vstack(expected).mean(0)
		bound = np.vstack(bound)
		active_tfs = np.vstack(active_tfs)

		# Calculate actual average promoter bound probabilities
		promoter_prob = (bound / np.vstack(promoters))
		promoter_prob[bound == 0] = 0
		actual = promoter_prob.mean(0)

		# Calculate how often the promoter binding was limited by the number of
		# active transcription factors
		active_limited_timesteps = (active_tfs == 0) & (bound > 0)
		frac_limited = active_limited_timesteps.sum(0) / active_limited_timesteps.shape[0]

		# Gene labels for limited points
		gene_ids = [tf_to_gene_id[tf] for tf in tf_ids]

		min_prob = min(expected[expected > 0].min(), actual[actual > 0].min())
		max_prob = np.max((expected, actual))

		# Create figure
		fig = plt.figure()
		ax = plt.gca()

		# Plot data
		ax.plot([0, max_prob], [0, max_prob], '--k')
		scatter = ax.scatter(expected, actual, c=frac_limited, alpha=0.5, cmap='hot', edgecolors='k')
		fig.colorbar(scatter)
		ax.set_xscale('symlog', linthreshx=min_prob)
		ax.set_yscale('symlog', linthreshy=min_prob)

		# Add text labels
		limited_cutoff = 0.05
		for tf, frac, x, y in zip(gene_ids, frac_limited, expected, actual):
			if frac < limited_cutoff:
				continue
			ax.text(x, 0, tf, ha='center', rotation=90, fontsize=6)
			ax.text(0, y, tf, va='center', fontsize=6)

		# Axes formatting
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)
		ax.set_xlabel(f'Expected promoter bound probability')
		ax.set_ylabel(f'Observed promoter bound probability')
		ax.set_title('Promoter binding probabilities for each transcription factor\n'
			'Color indicates fraction of time limited by the number of active TFs', fontsize=12)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
