"""Utilities for working with TOYA data"""
from __future__ import absolute_import, division, print_function

import re

import numpy as np
from typing import Any, Dict, Iterable, List, Tuple

from unum import Unum

from models.ecoli.processes.metabolism import (
	COUNTS_UNITS,
	VOLUME_UNITS,
	TIME_UNITS,
)
from six.moves import zip

FLUX_UNITS = COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS


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



def get_common_ids(toya_reaction_ids, root_to_id_indices_map):
	return list(sorted((
		set(toya_reaction_ids)
		& set(root_to_id_indices_map.keys())
	)))


def get_root_to_id_indices_map(sim_reaction_ids):
	# type: (np.ndarray) -> Dict[str, List[int]]
	"""Get a map from ID root to indices in sim_reaction_ids

	The reaction root is the portion of the ID that precedes the
	first underscore or space.

	Arguments:
		sim_reaction_ids: Simulation reaction IDs

	Returns:
		Map from ID root to a list of the indices in
		sim_reaction_ids of reactions having that root.
	"""

	raise NotImplementedError('This function does not provide all matches for certain reactions.'
		'Re-implement before using.')

	root_to_id_indices_map = dict()  # type: Dict[str, List[int]]
	matcher = re.compile("^([A-Za-z0-9-/.]+)")
	for i, rxn_id in enumerate(sim_reaction_ids):
		match = matcher.match(rxn_id)
		if match:
			root = match.group(1)
			root_to_id_indices_map.setdefault(root, []).append(i)
	return root_to_id_indices_map


def process_simulated_fluxes(
		output_ids, reaction_ids, reaction_fluxes):
	# type: (Iterable[str], Iterable[str], Unum) -> Tuple[Unum, Unum]
	"""Compute means and standard deviations of flux from simulation

	For a given output ID from output_ids, all reaction IDs from
	reaction_ids whose roots (before the first _ or space) match the output ID
	will have their data included in that output
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
	reaction_ids = np.array(reaction_ids)
	means = []  # type: List[np.ndarray]
	stdevs = []  # type: List[np.ndarray]
	for output_id in output_ids:
		time_course = []  # type: List[Unum]
		for i, rxn_id in enumerate(reaction_ids):
			if re.findall(output_id, rxn_id):
				reverse = -1 if re.findall("(reverse)", rxn_id) else 1
				matches = reaction_fluxes[:, i]
				if len(time_course):
					time_course += reverse * matches
				else:
					time_course = reverse * matches
		if len(time_course):
			time_course_ = time_course  # type: Any
			means.append(np.mean(time_course_).asNumber(FLUX_UNITS))
			stdevs.append(np.std(time_course_.asNumber(FLUX_UNITS)))
	means_ = FLUX_UNITS * np.array(means)
	stdevs_ = FLUX_UNITS * np.array(stdevs)
	return means_, stdevs_


def process_toya_data(output_ids, reaction_ids, data):
	# type: (Iterable[str], Iterable[str], Unum) -> Unum
	"""Filter toya fluxes or standard deviations by reaction ID

	Arguments:
		output_ids: IDs of reactions to include in the output.
		reaction_ids: IDs of the reactions whose fluxes and standard
			deviations are provided, in the order in which the
			rections' values appear in fluxes and stdevs.
		data: 1-dimensional numpy array (with units FLUX_UNITS)
			with average reaction fluxes or standard deviations.

	Returns:
		List of data in the order specified by output_ids. List is
		of units FLUX_UNITS.
	"""
	data_dict = dict(zip(reaction_ids, data))
	output_data = [
		data_dict[output_id].asNumber(FLUX_UNITS)
		for output_id in output_ids
	]
	return FLUX_UNITS * np.array(output_data)
