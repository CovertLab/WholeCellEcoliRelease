from __future__ import absolute_import, division, print_function

import os
from six.moves import cPickle
from typing import List, Optional, Tuple

import numpy as np
from matplotlib import pyplot as plt
from unum import Unum

from wholecell.io.tablereader import TableReader
from wholecell.utils import units, toya
from wholecell.analysis.plotting_tools import CMAP_COLORS_255
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot
from models.ecoli.processes.metabolism import (
	COUNTS_UNITS,
	VOLUME_UNITS,
	TIME_UNITS,
)

FLUX_UNITS = COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS

CMAP_COLORS = [[shade / 255. for shade in color] for color in CMAP_COLORS_255]
CMAP_OVER = [0, 1, 0.75]


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	"""Plot fluxome data"""

	@staticmethod
	def load_toya_data(validation_data_file, sim_data_file, sim_out_dir):
		# type: (str, str, str) -> Tuple[List[str], Unum, Unum]
		"""Load fluxome data from 2010 validation data

		Arguments:
			validation_data_file: Path to cPickle with validation data.
			sim_data_file: Path to cPickle with simulation data.
			sim_out_dir: Path to simulation output directory.

		Returns:
			Tuple of reaction IDs, fluxes, and standard deviations.
			Fluxes and standard deviations appear in the same order as
			their associated reactions do in the reaction ID list.
			Fluxes and standard deviations are numpy arrays with units
			FLUX_UNITS.
		"""
		validation_data = cPickle.load(open(validation_data_file, "rb"))
		sim_data = cPickle.load(open(sim_data_file, "rb"))
		cell_density = sim_data.constants.cell_density

		mass_listener = TableReader(os.path.join(sim_out_dir, "Mass"))
		cell_masses = mass_listener.readColumn("cellMass") * units.fg
		dry_masses = mass_listener.readColumn("dryMass") * units.fg
		mass_listener.close()

		toya_reactions = validation_data.reactionFlux.toya2010fluxes["reactionID"]
		toya_fluxes = toya.adjust_toya_data(
			validation_data.reactionFlux.toya2010fluxes["reactionFlux"],
			cell_masses,
			dry_masses,
			cell_density,
		)
		toya_stdevs = toya.adjust_toya_data(
			validation_data.reactionFlux.toya2010fluxes["reactionFluxStdev"],
			cell_masses,
			dry_masses,
			cell_density,
		)
		return toya_reactions, toya_fluxes, toya_stdevs

	@staticmethod
	def load_fba_data(simOutDir):
		# type: (str) -> Tuple[np.ndarray, Unum]
		"""Load fluxome balance analysis (FBA) data

		Arguments:
			simOutDir: Path to simulation output directory.

		Returns: Tuple of reaction IDs and, in the same order, the
			fluxes for each reaction. Fluxes array has units FLUX_UNITS.
		"""
		fba_results = TableReader(os.path.join(simOutDir, "FBAResults"))
		reaction_ids = np.array(fba_results.readAttribute("reactionIDs"))
		reaction_fluxes = FLUX_UNITS * np.array(
			fba_results.readColumn("reactionFluxes")
		)
		fba_results.close()
		return reaction_ids, reaction_fluxes

	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		# type: (str, str, str, str, str, Optional[dict]) -> None
		reaction_ids, reaction_fluxes = Plot.load_fba_data(simOutDir)
		toya_reactions, toya_fluxes, toya_stdevs = Plot.load_toya_data(
			validationDataFile, simDataFile, simOutDir)

		sim_flux_means, sim_flux_stdevs = toya.process_simulated_fluxes(
			toya_reactions, reaction_ids, reaction_fluxes
		)
		toya_flux_means = toya.process_toya_data(
			toya_reactions, toya_reactions, toya_fluxes)
		toya_flux_stdevs = toya.process_toya_data(
			toya_reactions, toya_reactions, toya_stdevs)

		correlation_coefficient = np.corrcoef(
			sim_flux_means.asNumber(FLUX_UNITS),
			toya_flux_means.asNumber(FLUX_UNITS),
		)[0, 1]
		plt.figure()

		plt.title("Central Carbon Metabolism Flux, Pearson R = {:.2}".format(
			correlation_coefficient))
		plt.errorbar(
			toya_flux_means.asNumber(FLUX_UNITS),
			sim_flux_means.asNumber(FLUX_UNITS),
			xerr=toya_flux_stdevs.asNumber(FLUX_UNITS),
			yerr=sim_flux_stdevs.asNumber(FLUX_UNITS),
			fmt="o", ecolor="k"
		)
		plt.ylabel("Mean WCM Reaction Flux {}".format(FLUX_UNITS.strUnit()))
		plt.xlabel("Toya 2010 Reaction Flux {}".format(FLUX_UNITS.strUnit()))

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
