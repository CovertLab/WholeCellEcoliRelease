from __future__ import absolute_import, division, print_function

import numpy as np

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)


def transcription_initiation_shuffle_params(sim_data, index):
	if index == 0:
		return CONTROL_OUTPUT, sim_data

	np.random.seed(index)
	idxs = np.arange(sim_data.process.transcription.rna_data.struct_array.shape[0])
	np.random.shuffle(idxs)
	sim_data.process.transcription.initiationShuffleIdxs = idxs

	return dict(
		shortName = "{}_transcriptionInitiationShuffle".format(index),
		desc = "Simulation of shuffled parameters for transcription initiation with seed {}.".format(index)
		), sim_data
