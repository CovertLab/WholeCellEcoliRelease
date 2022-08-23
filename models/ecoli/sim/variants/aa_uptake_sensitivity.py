"""
Add one amino acid to the minimal media condition and scale the uptake rate of
that amino acid but a range of factors.

Modifies:
	sim_data.external_state.saved_media
	sim_data.process.metabolism.max_specific_import_rates
	sim_data.process.metabolism.import_kcats_per_aa

Expected variant indices (dependent on length of FACTORS and sim_data.molecule_groups.amino_acids):
	0-7: add first amino acid to the media with different uptake rates
	8-167: add other amino acids with different uptake rates
"""

from .add_one_aa import add_one_aa


FACTORS = [0, 0.1, 0.5, 1, 1.5, 2, 5, 10]


def get_aa_index(index):
	return index // len(FACTORS)

def get_factor(index):
	return FACTORS[index % len(FACTORS)]

def aa_uptake_sensitivity(sim_data, index):
	# Use add_one_aa variant to add a specific amino acid to the media
	aa_idx = get_aa_index(index)
	_, sim_data = add_one_aa(sim_data, aa_idx)

	# Change the uptake rate for that amino acid to check the sensitivity
	factor = get_factor(index)
	sim_data.process.metabolism.max_specific_import_rates[aa_idx] *= factor
	sim_data.process.metabolism.import_kcats_per_aa[aa_idx] *= factor

	aa_id = sim_data.molecule_groups.amino_acids[aa_idx]

	return dict(
		shortName=f'{aa_id}:{factor}x',
		desc=f'Add {aa_id} to media at {factor}x uptake rate'
		), sim_data
