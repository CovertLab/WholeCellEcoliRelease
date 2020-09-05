"""
Variant to compare the impacts of reversing the orientation of rRNA genes.

Modifies:
	sim_data.process.transcription.rnaData["direction"]

Expected variant indices:
	0: control
	1: reverse orientation of all rRNA genes
"""

from __future__ import absolute_import, division, print_function

import numpy as np
from wholecell.utils import units

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)


def rrna_orientation(sim_data, index):
	if index == 1:
		is_rrna = sim_data.process.transcription.rna_data['is_rRNA']
		direction = sim_data.process.transcription.rna_data["direction"]
		coordinate = sim_data.process.transcription.rna_data[
			"replication_coordinate"]
		length = sim_data.process.transcription.rna_data["length"].asNumber(units.nt)

		# Reverse orientations of all rRNA genes
		coordinate[is_rrna] = coordinate[is_rrna] + np.multiply(
			length[is_rrna], direction[is_rrna])
		direction[is_rrna] = ~direction[is_rrna]

		return dict(
			shortName = "rrna_reversed",
			desc = "Simulation with all rRNA orientations reversed"
			), sim_data

	else:
		return CONTROL_OUTPUT, sim_data
