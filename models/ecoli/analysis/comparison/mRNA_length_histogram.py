"""
Compare histograms of mRNA lengths.
"""

from typing import Tuple

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import comparisonAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from validation.ecoli.validation_data import ValidationDataEcoli
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns
# noinspection PyUnresolvedReferences
from wholecell.io.tablereader import TableReader, TableReaderError
from wholecell.utils import units


FIGSIZE = (4, 4)
BOUNDS = [1, 5]
N_BINS = 20

class Plot(comparisonAnalysisPlot.ComparisonAnalysisPlot):
	def do_plot(self, reference_sim_dir, plotOutDir, plotOutFileName, input_sim_dir, unused, metadata):
		# noinspection PyUnusedLocal
		ap1, sim_data1, _ = self.setup(reference_sim_dir)
		# noinspection PyUnusedLocal
		ap2, sim_data2, _ = self.setup(input_sim_dir)

		def read_sims(ap, sim_data):
			rna_lengths = sim_data.process.transcription.rna_data['length'].asNumber(units.nt)
			is_mRNA = sim_data.process.transcription.rna_data['is_mRNA']

			mRNA_counts = read_stacked_columns(
				ap.get_cells(), 'mRNACounts', 'full_mRNA_counts')

			mean_mRNA_counts = mRNA_counts.mean(axis=0)
			mRNA_lengths = rna_lengths[is_mRNA]

			return mRNA_lengths, mean_mRNA_counts

		lengths1, mean_counts1 = read_sims(ap1, sim_data1)
		lengths2, mean_counts2 = read_sims(ap2, sim_data2)

		plt.figure(figsize=FIGSIZE)

		bins = np.logspace(BOUNDS[0], BOUNDS[1], N_BINS + 1)

		plt.hist(
			lengths1, bins=bins, weights=mean_counts1, alpha=0.5,
			label=f'reference ({np.mean(lengths1):.1f} $\pm$ {np.std(lengths1):.1f})')
		plt.hist(
			lengths2, bins=bins, weights=mean_counts2, alpha=0.5,
			label=f'input ({np.mean(lengths2):.1f} $\pm$ {np.std(lengths2):.1f})')

		plt.legend(prop={'size': 6})

		plt.xlim([10**1.5, 10**4.5])
		plt.xscale('log')
		plt.xlabel('mRNA lengths (nt)')
		plt.ylabel('mRNA counts')

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')

	def setup(self, inputDir: str) -> Tuple[
			AnalysisPaths, SimulationDataEcoli, ValidationDataEcoli]:
		"""Return objects used for analyzing multiple sims."""
		ap = AnalysisPaths(inputDir, variant_plot=True)
		sim_data = self.read_sim_data_file(inputDir)
		validation_data = self.read_validation_data_file(inputDir)
		return ap, sim_data, validation_data


if __name__ == "__main__":
	Plot().cli()
