import numpy as np
SCALING_FACTORS = np.arange(5, 16)

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def scalingFactorTotalIndices(sim_data):
	nScalingFactors = len(SCALING_FACTORS)
	nConditions = nScalingFactors + 1
	return nConditions


def scalingFactor(sim_data, index):
	nConditions = scalingFactorTotalIndices(sim_data)

	if index == 0:
		return CONTROL_OUTPUT, sim_data

	sim_data.scaling_factor = SCALING_FACTORS[index - 1]

	return dict(
		shortName = "{} factor".format(SCALING_FACTORS[index - 1]),
		desc = "Simulation uses factor of {}.".format(SCALING_FACTORS[index - 1])
		), sim_data
