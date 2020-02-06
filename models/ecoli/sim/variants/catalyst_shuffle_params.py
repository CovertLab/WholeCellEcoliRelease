CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)


def catalyst_shuffle_params(sim_data, index):
	if index == 0:
		return CONTROL_OUTPUT, sim_data

	import numpy as np
	np.random.seed(index)
	idxs = np.arange(len(sim_data.process.metabolism.reactions_with_catalyst))
	np.random.shuffle(idxs)
	sim_data.process.metabolism.catalystShuffleIdxs = idxs

	return dict(
		shortName = "{}_catalystShuffle".format(index),
		desc = "Simulation of shuffled parameters for catalysts in metabolism (to change zero/non-zero flux bounds) with seed {}.".format(index)
		), sim_data
