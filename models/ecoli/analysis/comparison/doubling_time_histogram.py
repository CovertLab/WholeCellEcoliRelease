"""
Compare histograms of doubling times.
"""

from typing import Tuple

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np
import os

from models.ecoli.analysis import comparisonAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from validation.ecoli.validation_data import ValidationDataEcoli
from wholecell.analysis.analysis_tools import exportFigure
# noinspection PyUnresolvedReferences
from wholecell.io.tablereader import TableReader, TableReaderError
from wholecell.utils import units


FIGSIZE = (4, 4)
DOUBLING_TIME_BOUNDS_MINUTES = [20, 180]
N_BINS = 32


class Plot(comparisonAnalysisPlot.ComparisonAnalysisPlot):
	def do_plot(self, reference_sim_dir, plotOutDir, plotOutFileName, input_sim_dir, unused, metadata):
		# noinspection PyUnusedLocal
		ap1, sim_data, _ = self.setup(reference_sim_dir)
		# noinspection PyUnusedLocal
		ap2, _, _ = self.setup(input_sim_dir)

		def read_sims(ap):
			sim_dirs = ap.get_cells()
			doubling_times_minutes = []

			for sim_dir in sim_dirs:
				try:
					sim_out_dir = os.path.join(sim_dir, "simOut")

					# Assume simulated time == doubling time
					time = TableReader(os.path.join(sim_out_dir, 'Main')).readColumn('time')
					doubling_time = time[-1] - time[0]
					doubling_times_minutes.append(doubling_time / 60.)

				except (TableReaderError, IndexError):
					continue

			return np.array(doubling_times_minutes)

		dt1 = read_sims(ap1)
		dt2 = read_sims(ap2)

		# Exclude sims that hit time limit
		dt1 = dt1[dt1 < 180]
		dt2 = dt2[dt2 < 180]

		if np.any(dt2 > 100):
			xlim = [20, 180]
		else:
			xlim = [20, 100]

		fig = plt.figure(figsize=FIGSIZE)
		ax = fig.add_subplot(1, 1, 1)

		bins = np.linspace(
			DOUBLING_TIME_BOUNDS_MINUTES[0],
			DOUBLING_TIME_BOUNDS_MINUTES[1],
			N_BINS + 1
			)

		ax.hist(
			dt1, bins=bins, alpha=0.5,
			label=f'reference (n={len(dt1)}, {np.mean(dt1):.1f} $\pm$ {np.std(dt1):.1f})')
		ax.hist(
			dt2, bins=bins, alpha=0.5,
			label=f'input (n={len(dt2)}, {np.mean(dt2):.1f} $\pm$ {np.std(dt2):.1f})')
		ax.legend(prop={'size': 8})

		ax.set_xlim(*xlim)

		ax.set_xlabel('Doubling time (min)')
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
