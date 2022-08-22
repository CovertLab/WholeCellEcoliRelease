"""
Compare histograms of mRNA masses.
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
BOUNDS = [0, 4]
N_BINS = 20

class Plot(comparisonAnalysisPlot.ComparisonAnalysisPlot):
	def do_plot(self, reference_sim_dir, plotOutDir, plotOutFileName, input_sim_dir, unused, metadata):
		# noinspection PyUnusedLocal
		ap1, sim_data1, _ = self.setup(reference_sim_dir)
		# noinspection PyUnusedLocal
		ap2, sim_data2, _ = self.setup(input_sim_dir)

		def read_sims(ap):
			# Ignore data from first two gens
			cell_paths = ap.get_cells(generation=np.arange(2, ap.n_generation))

			mRNA_mass = read_stacked_columns(
				cell_paths, 'Mass', 'mRnaMass',
				ignore_exception=True, fun=lambda x: x[0])

			return mRNA_mass

		mRNA_masses1 = read_sims(ap1)
		mRNA_masses2 = read_sims(ap2)

		fig = plt.figure(figsize=FIGSIZE)
		ax = fig.add_subplot(1, 1, 1)

		bins = np.linspace(BOUNDS[0], BOUNDS[1], N_BINS + 1)

		ax.hist(
			mRNA_masses1, bins=bins, alpha=0.5,
			label=f'reference ({np.mean(mRNA_masses1):.2f} $\pm$ {np.std(mRNA_masses1):.2f})')
		ax.axvline(np.mean(mRNA_masses1), ls='--', lw=2, c='C0')
		ax.hist(
			mRNA_masses2, bins=bins, alpha=0.5,
			label=f'input ({np.mean(mRNA_masses2):.2f} $\pm$ {np.std(mRNA_masses2):.2f})')
		ax.axvline(np.mean(mRNA_masses2), ls='--', lw=2, c='C1')

		ax.legend(prop={'size': 6})

		ax.set_xlim([0, 4])
		ax.set_xlabel('mRNA mass (fg)')
		ax.spines["top"].set_visible(False)
		ax.spines["right"].set_visible(False)

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
