from __future__ import absolute_import, division, print_function

import numpy as np

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)


def all_shuffle_params(sim_data, index):
	if index == 0:
		return CONTROL_OUTPUT, sim_data

	np.random.seed(index)

	# Shuffle transcription initiation
	idxs_transcriptionInitiation = np.arange(sim_data.process.transcription.rna_data.struct_array.shape[0])
	np.random.shuffle(idxs_transcriptionInitiation)
	sim_data.process.transcription.initiationShuffleIdxs = idxs_transcriptionInitiation

	# Shuffle translation efficiencies
	idxs_translationEfficiencies = np.arange(sim_data.process.translation.translation_efficiencies_by_monomer.shape[0])
	np.random.shuffle(idxs_translationEfficiencies)
	sim_data.process.translation.translationEfficienciesShuffleIdxs = idxs_translationEfficiencies

	# Shuffle monomer deg rates
	idxs_monomerDegRates = np.arange(len(sim_data.process.translation.monomer_data['deg_rate']))
	np.random.shuffle(idxs_monomerDegRates)
	sim_data.process.translation.monomerDegRateShuffleIdxs = idxs_monomerDegRates

	return dict(
		shortName = "{}_allShuffle".format(index),
		desc = "Simulation of shuffled parameters for all parameters implemented for shuffling (transcription initiation, translation efficiencies, monomer deg rates, metabolism kinetics shuffled with all reactions that have a catalyst) with seed {}.".format(index)
		), sim_data
