COEFF_VALUES = [0, 1e-9, 1e-6, 5e-5, 1e-5, 5e-4, 1e-4, 5e-3, 1e-3, 1e-2, 1e-1, .5]

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def metabolismKineticHomeostaticRatioTotalIndices(sim_data):
	nConditions = len(COEFF_VALUES) + 1
	return nConditions


def metabolismKineticHomeostaticRatio(sim_data, index):
	nConditions = metabolismKineticHomeostaticRatioTotalIndices(sim_data)

	if index % nConditions == 0:
		return CONTROL_OUTPUT, sim_data

	coeffValue = COEFF_VALUES[(index - 1) % nConditions]
	sim_data.constants.metabolismKineticObjectiveWeightLinear = coeffValue

	return dict(
		shortName = "kineticCoeff %.1g" % (coeffValue),
		desc = "Metabolic kinetic weighting coefficient set to %.3g." % (coeffValue)
		), sim_data