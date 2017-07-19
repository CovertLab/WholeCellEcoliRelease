
CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def translationEfficienciesShuffleParamsTotalIndices(sim_data):
	return None


def translationEfficienciesShuffleParams(sim_data, index):

	if index == 0:
		return CONTROL_OUTPUT, sim_data

	import numpy as np
	np.random.seed(index)
	idxs = np.arange(sim_data.process.translation.translationEfficienciesByMonomer.shape[0])
	np.random.shuffle(idxs)
	sim_data.process.translation.translationEfficienciesShuffleIdxs = idxs

	return dict(
		shortName = "{}_translationEfficienciesShuffle".format(index),
		desc = "Simulation of shuffled parameters for translation efficiencies with seed {}.".format(index)
		), sim_data
