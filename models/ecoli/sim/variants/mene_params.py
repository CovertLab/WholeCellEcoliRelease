"""
Variant to compare varying levels of menE expression to show subgenerational
expression.

Modifies:
	sim_data.process.transcription_regulation.basal_prob
	sim_data.process.transcription_regulation.delta_prob
	sim_data.process.transcription.rnaExpression

Expected variant indices (dependent on FACTORS):
	4: control
	0 (low expression) - 8 (high expression)
"""

from __future__ import absolute_import, division, print_function

import numpy as np

FACTORS = [0.1, 0.125, 0.25, 0.5, 1, 2., 4., 8., 10.]

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)


def mene_params(sim_data, index):
	if index == FACTORS.index(1):
		return CONTROL_OUTPUT, sim_data

	adjustFactor = FACTORS[index]
	mene_rna_index = np.where(sim_data.process.transcription.rnaData["id"] == "EG12437_RNA[c]")[0]
	mene_delta_prob_indices = np.where(sim_data.process.transcription_regulation.delta_prob["deltaI"] == mene_rna_index)[0]

	# Adjust transcript synthesis probability in basal_prob and delta_prob
	sim_data.process.transcription_regulation.basal_prob[mene_rna_index] *= adjustFactor
	sim_data.process.transcription_regulation.delta_prob["deltaV"][mene_delta_prob_indices] *= adjustFactor

	# Adjust rna expression data in order to impact initial condition
	sim_data.process.transcription.rnaExpression[sim_data.condition][mene_rna_index] *= adjustFactor

	return dict(
		shortName = "{}_meneParams".format(index),
		desc = "Simulation with menE synthesis probability increased by the factor {}.".format(FACTORS[index])
		), sim_data
