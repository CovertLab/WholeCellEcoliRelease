"""
Variant to compare varying levels of gene expression to show subgenerational
expression.  Useful with subgen_expression variant analysis.

Modifies:
	sim_data.process.transcription_regulation.recruitmentData
	sim_data.process.transcription.rnaExpression
	sim_data.process.transcription.rnaSynthProb

Expected variant indices (dependent on FACTORS):
	4: control
	0 (low expression) - 8 (high expression)
"""

import numpy as np


FACTORS = [0.1, 0.125, 0.25, 0.5, 1, 2., 4., 8., 10.]
RNA_ID = 'EG10683_RNA[c]'  # pabB
CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)


def subgen_expression_total_indices(sim_data):
	return len(FACTORS)

def subgen_expression(sim_data, index):
	if index == FACTORS.index(1):
		return CONTROL_OUTPUT, sim_data

	factor = FACTORS[index]
	geneIndex = np.where(sim_data.process.transcription.rnaData["id"] == RNA_ID)[0][0]

	recruitment_mask = np.array([hi == geneIndex
		for hi in sim_data.process.transcription_regulation.recruitmentData['hI']])
	for synth_prob in sim_data.process.transcription.rnaSynthProb.values():
		synth_prob[geneIndex] *= factor
	for exp in sim_data.process.transcription.rnaExpression.values():
		exp[geneIndex] *= factor
	sim_data.process.transcription_regulation.recruitmentData['hV'][recruitment_mask] *= factor

	# Renormalize parameters
	for synth_prob in sim_data.process.transcription.rnaSynthProb.values():
		synth_prob /= synth_prob.sum()
	for exp in sim_data.process.transcription.rnaExpression.values():
		exp /= exp.sum()

	return dict(
		shortName = "{}_subgen_expression".format(index),
		desc = "Simulation with pabB synthesis probability adjusted by the factor {}.".format(FACTORS[index])
		), sim_data
