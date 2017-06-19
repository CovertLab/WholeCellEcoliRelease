
CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def transcriptionInitiationShuffleParamsTotalIndices(sim_data):
	return None


def transcriptionInitiationShuffleParams(sim_data, index):

	if index == 0:
		return CONTROL_OUTPUT, sim_data

	import numpy as np
	np.random.seed(index)
	idxs = np.arange(sim_data.process.transcription.rnaData.struct_array.shape[0])
	np.random.shuffle(idxs)
	sim_data.process.transcription.initiationShuffleIdxs = idxs

	return dict(
		shortName = "{}_transcriptionInitiationShuffle".format(index),
		desc = "Simulation of shuffled parameters for transcription initiation with seed {}.".format(index)
		), sim_data
