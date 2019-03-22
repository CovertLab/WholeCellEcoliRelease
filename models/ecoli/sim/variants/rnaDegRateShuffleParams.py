
CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def rnaDegRateShuffleParamsTotalIndices(sim_data):
	return None


def rnaDegRateShuffleParams(sim_data, index):

	if index == 0:
		return CONTROL_OUTPUT, sim_data

	import numpy as np
	np.random.seed(index)
	idxs = np.arange(len(sim_data.process.transcription.rnaData["degRate"]))
	np.random.shuffle(idxs)
	sim_data.process.transcription.rnaDegRateShuffleIdxs = idxs

	return dict(
		shortName = "{}_monomerDegRateShuffle".format(index),
		desc = "Simulation of shuffled parameters for monomer degradation rates with seed {}.".format(index)
		), sim_data
