import numpy as np

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)


def kinetic_catalyst_shuffle_params(sim_data, index):
	if index == 0:
		return CONTROL_OUTPUT, sim_data

	np.random.seed(index)
	idxs = np.arange(len(sim_data.process.metabolism.reactions_with_catalyst))
	np.random.shuffle(idxs)
	nTargets = len(sim_data.process.metabolism.kinetic_constraint_reactions)
	sim_data.process.metabolism.kineticTargetShuffleRxns = [sim_data.process.metabolism.reactions_with_catalyst[idx] for idx in idxs[:nTargets]]

	return dict(
		shortName = "{}_kineticTargetShuffle".format(index),
		desc = "Simulation of shuffled parameters for metabolism kinetics with all reactions with a catalyst with seed {}.".format(index)
		), sim_data
