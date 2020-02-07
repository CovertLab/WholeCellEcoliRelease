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
	idxs_transcriptionInitiation = np.arange(sim_data.process.transcription.rnaData.struct_array.shape[0])
	np.random.shuffle(idxs_transcriptionInitiation)
	sim_data.process.transcription.initiationShuffleIdxs = idxs_transcriptionInitiation

	# Shuffle translation efficiencies
	idxs_translationEfficiencies = np.arange(sim_data.process.translation.translationEfficienciesByMonomer.shape[0])
	np.random.shuffle(idxs_translationEfficiencies)
	sim_data.process.translation.translationEfficienciesShuffleIdxs = idxs_translationEfficiencies

	# Shuffle monomer deg rates
	idxs_monomerDegRates = np.arange(len(sim_data.process.translation.monomerData["degRate"]))
	np.random.shuffle(idxs_monomerDegRates)
	sim_data.process.translation.monomerDegRateShuffleIdxs = idxs_monomerDegRates

	# Shuffle kinetic catalysts
	idxs_kineticCatalysts = np.arange(len(sim_data.process.metabolism.reactions_with_catalyst))
	np.random.shuffle(idxs_kineticCatalysts)
	nTargets = len(sim_data.process.metabolism.kinetic_constraint_reactions)
	sim_data.process.metabolism.kineticTargetShuffleRxns = [sim_data.process.metabolism.reactions_with_catalyst[idx] for idx in idxs_kineticCatalysts[:nTargets]]

	return dict(
		shortName = "{}_allShuffle".format(index),
		desc = "Simulation of shuffled parameters for all parameters implemented for shuffling (transcription initiation, translation efficiencies, monomer deg rates, metabolism kinetics shuffled with all reactions that have a catalyst) with seed {}.".format(index)
		), sim_data
