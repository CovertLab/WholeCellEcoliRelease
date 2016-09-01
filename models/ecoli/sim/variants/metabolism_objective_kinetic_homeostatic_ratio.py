COEFF_VALUES = [0, 1e-9, 1e-6, 1e-3, 1e-2, 1e-1, .25, .5, .75, (1-1e-1), (1-1e-2), (1-1e-3), (1-1e-6), (1-1e-9), 1]

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
	sim_data.constants.metabolismKineticObjectiveWeight = coeffValue

	return dict(
		shortName = "kineticCoeff %.1f" % (coeffValue),
		desc = "Metabolic kinetic weighting coefficient set to %.2f." % (coeffValue)
		), sim_data