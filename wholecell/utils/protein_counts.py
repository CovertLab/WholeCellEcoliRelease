"""Utilities for counting simulated proteins"""


from __future__ import absolute_import, division

import numpy as np


def get_sim_wisniewski_counts(
	monomer_counts, wisniewski_monomer_ids, sim_monomer_ids
):
	"""Get simulated monomer counts that correspond to Wisniewski data

	Arguments:
		monomer_counts: Simulated monomer counts (from translation
			process).
		wisniewski_monomer_ids: Monomer IDs from Wisniewski data. IDs
			must appear in same order as in Wisniewski monomer counts.
		sim_monomer_ids: IDs of monomers in the same order as
			monomer_counts.

	Returns:
		The simulated counts of the monomers that appear in the
		Wisniewski data, in the same order as in the Wisniewski data.
	"""
	# type: (np.ndarray, np.ndarray, np.ndarray) -> np.ndarray
	sim_ids_lst = sim_monomer_ids.tolist()
	wisniewski_ids_lst = wisniewski_monomer_ids.tolist()
	wisniewski_idx = [sim_ids_lst.index(x) for x in wisniewski_ids_lst]
	avg_counts = monomer_counts.mean(axis=0)
	return avg_counts[wisniewski_idx]


def get_sim_schmidt_counts(
	monomer_counts, schmidt_monomer_ids, sim_monomer_ids
):
	"""Get simulated monomer counts that correspond to Schmidt data

	Arguments:
		monomer_counts: Simulated monomer counts (from translation
			process).
		schmidt_monomer_ids: Monomer IDs from Schmidt data. IDs
			must appear in same order as in Schmidt monomer counts.
		sim_monomer_ids: IDs of monomers in the same order as
			monomer_counts.

	Returns:
		The simulated counts of the monomers that appear in the
		Schmidt data, in the same order as in the Schmidt data.
	"""
	# type: (np.ndarray, np.ndarray, np.ndarray) -> np.ndarray
	sim_ids_lst = sim_monomer_ids.tolist()
	schmidt_ids_lst = schmidt_monomer_ids.tolist()
	schmidt_idx = [sim_ids_lst.index(x) for x in schmidt_ids_lst]
	avg_counts = monomer_counts.mean(axis=0)
	return avg_counts[schmidt_idx]
