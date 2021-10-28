"""
Vary amino acid synthesis network parameters to see the effects on growth rate
and elongation rate in different media conditions.

Modifies:
	Attributes from add_one_aa variant
	sim_data.process.metabolism.aa_kcats
	sim_data.process.metabolism.aa_kis
	sim_data.process.metabolism.aa_upstream_kms
	sim_data.process.metabolism.aa_reverse_kms
	sim_data.process.metabolism.aa_degradation_kms

Expected variant indices (dependent on length of FACTORS, PARAMETERS,
sim_data.molecule_groups.amino_acids and MEDIA_IDS as noted below each group):
	0-5: range of values for first parameter, first amino acid, first media condition
		# FACTORS
	0-35: range of parameters and values for first amino acid, first media condition
		# FACTORS x # PARAMETERS
	0-755: range of parameters and values over all amino acids for first media condition
		# FACTORS x # PARAMETERS x # amino acids
	0-1511: all changes
		# FACTORS x # PARAMETERS x # amino acids x # MEDIA_IDS

TODO:
	- Run this for all AA additions or just Glt and control?
"""

from .add_one_aa import add_one_aa

import numpy as np


FACTORS = [0, 0.1, 0.5, 2, 5, 10]  # TODO: run factor of 1 once for each media condition
PARAMETERS = ['aa_kcats_fwd', 'aa_kcats_rev', 'aa_kis', 'aa_upstream_kms', 'aa_reverse_kms', 'aa_degradation_kms']
N_PARAM_VALUES = len(FACTORS) * len(PARAMETERS)
MEDIA_IDS = [5, 19]  # Glt and control for now (need to update the analysis plot if changed)


def get_n_aas(sim_data):
	"""Get the number of amino acids to vary parameters for in case it is not all of them"""
	return len(sim_data.molecule_groups.amino_acids)

def get_media_index(index, sim_data):
	n_aas = get_n_aas(sim_data)
	return MEDIA_IDS[index // (N_PARAM_VALUES * n_aas)]

def get_aa_index(index, sim_data):
	n_aas = get_n_aas(sim_data)
	sub_index = index % (N_PARAM_VALUES * n_aas)
	return sub_index // N_PARAM_VALUES

def get_adjustment(index):
	sub_index = index % N_PARAM_VALUES
	param_index = sub_index // len(FACTORS)
	factor_index = sub_index % len(FACTORS)
	return PARAMETERS[param_index], FACTORS[factor_index]

def aa_synthesis_sensitivity(sim_data, index):
	# Use add_one_aa variant to add a specific amino acid to the media
	media_idx = get_media_index(index, sim_data)
	_, sim_data = add_one_aa(sim_data, media_idx)

	# Change the uptake rate for that amino acid to check the sensitivity
	aa_idx = get_aa_index(index, sim_data)
	param, factor = get_adjustment(index)
	values = getattr(sim_data.process.metabolism, param)
	if np.all(~np.isfinite(values[aa_idx])) or np.all(values[aa_idx] == 0):
		# Skip sims for parameters that are 0 or inf and will not be updated
		raise ValueError('No change to params - not running variant sims.')
	values[aa_idx] *= factor

	aa_id = sim_data.molecule_groups.amino_acids[aa_idx]
	media_id = sim_data.molecule_groups.amino_acids[media_idx]

	return dict(
		shortName=f'{media_id}:{aa_id} {param} {factor}x',
		desc=f'{media_id} added to media: adjusted {aa_id} {param} by {factor}x'
		), sim_data
