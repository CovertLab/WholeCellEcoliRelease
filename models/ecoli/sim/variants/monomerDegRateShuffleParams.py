
CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def monomerDegRateShuffleParamsTotalIndices(sim_data):
	return None


def monomerDegRateShuffleParams(sim_data, index):

	if index == 0:
		return CONTROL_OUTPUT, sim_data

	import numpy as np
	np.random.seed(index)
	idxs = np.arange(len(sim_data.process.translation.monomerData["degRate"]))
	np.random.shuffle(idxs)
	sim_data.process.translation.monomerDegRateShuffleIdxs = idxs

	return dict(
		shortName = "{}_monomerDegRateShuffle".format(index),
		desc = "Simulation of shuffled parameters for monomer degradation rates with seed {}.".format(index)
		), sim_data
