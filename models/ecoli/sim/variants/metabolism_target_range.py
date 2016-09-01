COEFF_VALUES = [-2,-1,-.5, -.1, -.01, .01, .1, .5, 1, 2]

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def metabolismTargetRangeTotalIndices(sim_data):
	nGenes = sim_data.process.transcription.rnaData.fullArray().size
	nConditions = len(COEFF_VALUES) + 1
	return nConditions


def metabolismTargetRange(sim_data, index):
	nConditions = metabolismTargetRangeTotalIndices(sim_data)

	if index % nConditions == 0:
		return CONTROL_OUTPUT, sim_data

	coeffValue = COEFF_VALUES[(index - 1) % nConditions]
	sim_data.constants.metabolismTargetRangeConstant = coeffValue

	return dict(
		shortName = "rangeCoeff %.1f" % (coeffValue),
		desc = "Metabolic target range coefficient set to %.2f." % (coeffValue)
		), sim_data