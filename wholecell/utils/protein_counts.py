"""Utilities for counting simulated proteins"""

from __future__ import absolute_import, division, print_function

from typing import cast, List, Tuple

import numpy as np


def get_simulated_validation_counts(
		validation_counts, monomer_counts, validation_ids, simulation_ids):
	# type: (np.ndarray, np.ndarray, np.ndarray, np.ndarray) -> Tuple[np.ndarray, np.ndarray]
	"""
	Get simulated counts and validation counts of monomers that exist in both
	the simulation and the validation data

	Arguments:
		validation_counts: Monomer counts from validation data.
		monomer_counts: Simulated monomer counts (from translation
			process).
		validation_ids: Monomer IDs from validation data. IDs
			must appear in same order as in validation_counts.
		simulation_ids: IDs of monomers in the same order as
			monomer_counts.

	Returns:
		The simulated counts of the monomers that appear in the
		validation data, and the validation counts of the monomers in the same
		order.
	"""
	avg_sim_counts = monomer_counts.mean(axis=0)

	sim_ids_lst = cast(List[str], simulation_ids.tolist())
	val_ids_lst = cast(List[str], validation_ids.tolist())
	overlapping_ids_set = set(sim_ids_lst) & set(val_ids_lst)

	sim_id_to_index_map = {
		sim_id: i for i, sim_id in enumerate(sim_ids_lst)
		if sim_id in overlapping_ids_set}
	val_id_to_index_map = {
		val_id: i for i, val_id in enumerate(val_ids_lst)
		if val_id in overlapping_ids_set}

	overlapping_ids_list = list(overlapping_ids_set)
	sim_filtered_idx = np.array([
		sim_id_to_index_map[monomer_id] for monomer_id in overlapping_ids_list
		])
	val_filtered_idx = np.array([
		val_id_to_index_map[monomer_id] for monomer_id in overlapping_ids_list
		])

	return avg_sim_counts[sim_filtered_idx], validation_counts[val_filtered_idx]
