
CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def kineticCatalystShuffleParamsTotalIndices(sim_data):
	return None


def kineticCatalystShuffleParams(sim_data, index):

	if index == 0:
		return CONTROL_OUTPUT, sim_data

	import numpy as np
	np.random.seed(index)
	idxs = np.arange(len(sim_data.process.metabolism.reactionCatalystsList))
	np.random.shuffle(idxs)
	nTargets = len(sim_data.process.metabolism.constrainedReactionList)
	sim_data.process.metabolism.kineticTargetShuffleRxns = [sim_data.process.metabolism.reactionCatalystsList[idx] for idx in idxs[:nTargets]]

	return dict(
		shortName = "{}_kineticTargetShuffle".format(index),
		desc = "Simulation of shuffled parameters for metabolism kinetics with all reactions with a catalyst with seed {}.".format(index)
		), sim_data
