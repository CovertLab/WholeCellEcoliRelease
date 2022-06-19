"""
Shows dynamics of attenuated genes and probabilities associated with attenuation.
Expected and actual probabilities should converge by the end of a generation
if the probability is mostly constant throughout the sim.  Useful to see which
RNA are being attenuated and if that causes a reduction in the number of RNAs
transcribed.

TODO:
- add totals for all attenuated genes
- cleanup clutter (remove duplicate axis labels, legends, titles)
- gene symbol for RNA IDs
"""

import os

from matplotlib import pyplot as plt
from matplotlib import gridspec
import numpy as np

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	_suppress_numpy_warnings = True

	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Listeners used
		transcription_reader = TableReader(os.path.join(simOutDir, 'TranscriptElongationListener'))

		# Load data
		rna_ids = transcription_reader.readAttribute('rnaIds')
		attenuated_rnas = transcription_reader.readAttribute('attenuated_rnas')
		n_attenuated = len(attenuated_rnas)
		attenuated_mask = np.array([rna in set(attenuated_rnas) for rna in rna_ids])

		sim_time = transcription_reader.readColumn('time')
		expected_probability = transcription_reader.readColumn('attenuation_probability')
		counts_attenuated = transcription_reader.readColumn('counts_attenuated')
		counts_synthesized = transcription_reader.readColumn('countRnaSynthesized')[:, attenuated_mask]

		cum_attenuated = np.cumsum(counts_attenuated, axis=0)
		cum_synthesized = np.cumsum(counts_synthesized, axis=0)
		actual_probability = cum_attenuated / (cum_attenuated + cum_synthesized)

		plt.figure(figsize=(20, 20))
		n_cols = 4
		n_subcols = 2
		gs = gridspec.GridSpec(nrows=int(np.ceil(n_attenuated/n_cols)), ncols=n_subcols*n_cols)

		for i in range(n_attenuated):
			row = i // n_cols
			col = i % n_cols

			ax = plt.subplot(gs[row, n_subcols*col])
			ax.plot(sim_time, cum_attenuated[:, i], label='Attenuated')
			ax.plot(sim_time, cum_synthesized[:, i], label='Synthesized')
			ax.set_ylabel('Cumulative counts', fontsize=6)
			ax.tick_params(labelsize=6)
			ax.set_title(attenuated_rnas[i], fontsize=6)
			plt.legend(fontsize=6)

			ax = plt.subplot(gs[row, n_subcols*col+1])
			ax.plot(sim_time, expected_probability[:, i], label='Expected')
			ax.plot(sim_time, actual_probability[:, i], label='Actual')
			ax.set_ylabel('Probability', fontsize=6)
			ax.tick_params(labelsize=6)
			ax.set_title(attenuated_rnas[i], fontsize=6)
			plt.legend(fontsize=6)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
