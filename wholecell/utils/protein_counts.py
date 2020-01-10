"""Utilities for counting simulated proteins"""


from __future__ import absolute_import, division

import numpy as np


def get_simulated_validation_counts(
	monomer_counts, validation_ids, simulation_ids
):
	"""Get simulated monomer counts that correspond to Wisniewski data

	Arguments:
		monomer_counts: Simulated monomer counts (from translation
			process).
		validation_ids: Monomer IDs from validation data. IDs
			must appear in same order as in validation monomer counts.
		simulation_ids: IDs of monomers in the same order as
			monomer_counts.

	Returns:
		The simulated counts of the monomers that appear in the
		validation data, in the same order as in the validation data.
	"""
	# type: (np.ndarray, np.ndarray, np.ndarray) -> np.ndarray
	sim_ids_lst = simulation_ids.tolist()
	validation_ids_lst = validation_ids.tolist()
	validation_idx = [sim_ids_lst.index(x) for x in validation_ids_lst]
	avg_counts = monomer_counts.mean(axis=0)
	return avg_counts[validation_idx]
