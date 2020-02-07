import numpy as np

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)


def kinetic_target_shuffle_params(sim_data, index):
	if index == 0:
		return CONTROL_OUTPUT, sim_data

	np.random.seed(index)
	idxs = np.arange(len(sim_data.process.metabolism.kinetic_constraint_reactions))
	np.random.shuffle(idxs)
	sim_data.process.metabolism.kineticTargetShuffleIdxs = idxs

	return dict(
		shortName = "{}_kineticTargetShuffle".format(index),
		desc = "Simulation of shuffled parameters for metabolism kinetics with seed {}.".format(index)
		), sim_data
