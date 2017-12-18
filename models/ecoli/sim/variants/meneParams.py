import numpy as np

FACTORS = [0.1, 0.125, 0.25, 0.5, 1, 2., 4., 8., 10.]

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def meneParamsTotalIndices(sim_data):
	return None

def meneParams(sim_data, index):
	if index == FACTORS.index(1):
		return CONTROL_OUTPUT, sim_data

	adjustFactor = FACTORS[index]

	# Adjust transcript synthesis probability in recruitmentData
	mene_rnaIndex = np.where(sim_data.process.transcription.rnaData["id"] == "EG12437_RNA[c]")[0]
	mene_recruitmentIndices = np.where(sim_data.process.transcription_regulation.recruitmentData["hI"] == mene_rnaIndex)[0]
	sim_data.process.transcription_regulation.recruitmentData["hV"][mene_recruitmentIndices] *= adjustFactor

	# Adjust rna expression data in order to impact initial condition
	sim_data.process.transcription.rnaExpression[sim_data.condition][mene_rnaIndex] *= adjustFactor

	return dict(
		shortName = "{}_meneParams".format(index),
		desc = "Simulation with menE synthesis probability increased by the factor {}.".format(FACTORS[index])
		), sim_data
