from __future__ import absolute_import, division, print_function

import numpy as np

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)


def rna_deg_rate_shuffle_params(sim_data, index):
	if index == 0:
		return CONTROL_OUTPUT, sim_data

	np.random.seed(index)
	idxs = np.arange(len(sim_data.process.transcription.rna_data['deg_rate']))
	np.random.shuffle(idxs)
	sim_data.process.transcription.rnaDegRateShuffleIdxs = idxs

	return dict(
		shortName = "{}_monomerDegRateShuffle".format(index),
		desc = "Simulation of shuffled parameters for monomer degradation rates with seed {}.".format(index)
		), sim_data
