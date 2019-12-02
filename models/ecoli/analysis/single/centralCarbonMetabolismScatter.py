"""
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/29/2016
"""

from __future__ import absolute_import
from __future__ import division

import os
import cPickle
import re
from typing import List, Iterable, Tuple

import numpy as np
from matplotlib import pyplot as plt
from unum import Unum

from wholecell.io.tablereader import TableReader
from wholecell.utils import units
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
			validationDataFile: Path to cPickle with validation data.
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
		cell_density = sim_data.constants.cellDensity

		mass_listener = TableReader(os.path.join(sim_out_dir, "Mass"))
		cell_masses = mass_listener.readColumn("cellMass") * units.fg
		dry_masses = mass_listener.readColumn("dryMass") * units.fg
		mass_listener.close()

		toya_reactions = validation_data.reactionFlux.toya2010fluxes["reactionID"]
		toya_fluxes = Plot.adjust_toya_data(
			validation_data.reactionFlux.toya2010fluxes["reactionFlux"],
			cell_masses,
			dry_masses,
			cell_density,
		)
		toya_stdevs = Plot.adjust_toya_data(
			validation_data.reactionFlux.toya2010fluxes["reactionFluxStdev"],
			cell_masses,
			dry_masses,
			cell_density,
		)
		return toya_reactions, toya_fluxes, toya_stdevs

	@staticmethod
	def adjust_toya_data(data, cell_masses, dry_masses, cell_density):
		# type: (np.ndarray, np.ndarray, np.ndarray, float) -> Unum
		"""Adjust fluxes or stdevs for dry mass fraction, density, units.

		The provided data are multiplied by the average dry mass
		fraction and by the cell density.

		Arguments:
			data: Flux or standard deviation data vector to process.
			cell_masses: Vector of cell mass over time. Units should be
				stored elementwise.
			dry_masses: Vector of cell dry mass over time. Units should
				be stored elementwise.
			cell_density: Constant density of cell.

		Returns:
			Vector of adjusted data. Units are applied to the whole
			vector.
		"""
		dry_mass_frac_average = np.mean(dry_masses / cell_masses)
		return FLUX_UNITS * np.array([
			(dry_mass_frac_average * cell_density * x).asNumber(FLUX_UNITS)
			for x in data
		])

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

	def do_plot(
		self,
		simOutDir,  # type: str
		plotOutDir,  # type: str
		plotOutFileName,  # type: str
		simDataFile,  # type: str
		validationDataFile,  # type: str
		metadata
	):
		if not os.path.isdir(simOutDir):
			raise Exception(
				"simOutDir does not currently exist as a directory")

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		reaction_ids, reaction_fluxes = Plot.load_fba_data(simOutDir)
		toya_reactions, toya_fluxes, toya_stdevs = Plot.load_toya_data(
			validationDataFile, simDataFile, simOutDir)
		common_ids = [
			rxn_id for rxn_id in toya_reactions
			if Plot.regex_in_list(rxn_id, reaction_ids)
		]

		sim_flux_means, sim_flux_stdevs = Plot.process_simulated_fluxes(
			common_ids, reaction_ids, reaction_fluxes)
		toya_flux_means, toya_flux_stdevs = Plot.process_toya_fluxes(
			common_ids, toya_reactions, toya_fluxes, toya_stdevs)

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

	@staticmethod
	def regex_in_list(regex, lst):
		# type: (str, List[str]) -> bool
		"""Search a list for a regex match.

		Arguments:
			regex: Regular expression to search for matches for.
			lst: List to search.

		Returns:
			True if and only if a match is found.
		"""
		matcher = re.compile(regex)
		for elem in lst:
			if matcher.match(elem):
				return True
		return False

	@staticmethod
	def process_simulated_fluxes(
		output_ids,  # type: Iterable[str]
		reaction_ids,  # type: Iterable[str]
		reaction_fluxes,  # type: Unum
	):
		# type: (...) -> Tuple[Unum, Unum]
		"""Compute means and standard deviations of flux from simulation

		For a given output ID from output_ids, all reaction IDs from
		reaction_ids that match the output ID (treating output ID as a
		regular expression) will have their data included in that output
		ID's mean and standard deviation.

		Arguments:
			output_ids: IDs of reactions to include in output
			reaction_ids: IDs of the reactions in the order that they
				appear in reaction_fluxes
			reaction_fluxes: 2-dimensional matrix of reaction fluxes
				where each column corresponds to a reaction (in the
				order specified by reaction_ids) and each row is a time
				point. Should have units FLUX_UNITS and be a numpy
				matrix.

		Returns:
			Tuple of the lists of mean fluxes and standard deviations for each
			reaction ID in output_ids. List elements are in the same
			order as their associated reaction IDs in output_ids. Both
			lists will have units FLUX_UNITS.
		"""
		means = []
		stdevs = []
		for output_id in output_ids:
			time_course = []
			for rxn_id in reaction_ids:
				if not re.findall(output_id, rxn_id):
					continue
				reverse = -1 if re.findall("(reverse)", rxn_id) else 1
				matches = reaction_fluxes[:, np.where(reaction_ids == rxn_id)]
				if time_course:
					time_course += reverse * matches
				else:
					time_course = reverse * matches
			if time_course:
				means.append(np.mean(time_course).asNumber(FLUX_UNITS))
				stdevs.append(np.std(time_course.asNumber(FLUX_UNITS)))
		means = FLUX_UNITS * np.array(means)
		stdevs = FLUX_UNITS * np.array(stdevs)
		return means, stdevs

	@staticmethod
	def process_toya_fluxes(
		output_ids,  # type: Iterable[str]
		reaction_ids,  # type: Iterable[str]
		fluxes,  # type: Unum
		stdevs,  # type: Unum
	):
		# type: (...) -> Tuple[Unum, Unum]
		"""Filter toya fluxes and standard deviations by reaction ID

		Arguments:
			output_ids: IDs of reactions to include in the output.
			reaction_ids: IDs of the reactions whose fluxes and standard
				deviations are provided, in the order in which the
				rections' values appear in fluxes and stdevs.
			fluxes: 1-dimensional numpy array (with units FLUX_UNITS)
				with average reaction fluxes.
			stdevs: 1-dimensional numpy array (with units FLUX_UNITS)
				with the standard deviations of reaction fluxes.

		Returns:
			Tuple of the reaction fluxes and standard deviations, with
			each list in the order specified by output_ids. Both lists
			are of units FLUX_UNITS.
		"""
		fluxes_dict = dict(zip(reaction_ids, fluxes))
		stdevs_dict = dict(zip(reaction_ids, stdevs))
		output_fluxes = [
			fluxes_dict[output_id].asNumber(FLUX_UNITS)
			for output_id in output_ids
		]
		output_stdevs = [
			stdevs_dict[output_id].asNumber(FLUX_UNITS)
			for output_id in output_ids
		]
		output_fluxes = FLUX_UNITS * np.array(output_fluxes)
		output_stdevs = FLUX_UNITS * np.array(output_stdevs)
		return output_fluxes, output_stdevs


if __name__ == "__main__":
	Plot().cli()
