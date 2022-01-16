"""
Remove amino acids from the minimal media plus amino acids condition after 10 minutes.

Modifies:
	sim_data.condition
	sim_data.external_state.saved_media
	sim_data.external_state.current_timeline_id
	sim_data.external_state.saved_timelines
	sim_data.expectedDryMassIncreaseDict
	sim_data.nutrient_to_doubling_time
	sim_data.translation_supply_rate
	sim_data.process.translation.ribosomeElongationRateDict
	sim_data.process.translation.ribosomeFractionActiveDict
	sim_data.process.transcription.rnaPolymeraseElongationRateDict
	sim_data.process.transcription.rnaSynthProbFraction
	sim_data.process.transcription.rnaSynthProbRProtein
	sim_data.process.transcription.rnaSynthProbRnaPolymerase
	sim_data.process.transcription.rnapFractionActiveDict

Expected variant indices:
	0: All amino acids
	1: Shift to 12 amino acids
	2: Shift to 6 amino acids
	3: Minimal media (no amino acids)
	4-22: Shift to single amino acids
	23: Shift to no amino acids
"""

import numpy as np


SHIFT_TIME = 600  # 10 minutes

# For comparison to Taheri-Araghi et al. Cell-Size Control and Homeostasis in Bacteria. 2015.
AA_6 = ['MET', 'HIS', 'ARG', 'PRO', 'THR', 'TRP']
AA_12 = AA_6 + ['SER', 'LEU', 'TYR', 'L-ALPHA-ALANINE', 'ASN', 'L-ASPARTATE']
N_COMBOS = 2  # Number of amino acid combo variants (eg AA_6, AA_12)
N_FIXED_VARIANTS = N_COMBOS + 2  # +2 for rich and minimal variants

SKIP_AAS = {
	'L-SELENOCYSTEINE',  # Required to be in media
	'CYS',  # Import not currently implemented
	}


def remove_aas_shift(sim_data, index):
	# Remove amino acids from the media
	condition_label = 'with_aa'
	base_media_id = 'minimal_plus_amino_acids'
	minimal_media_id = 'minimal'
	aa_ids = [aa[:-3] for aa in sim_data.molecule_groups.amino_acids if aa[:-3] not in SKIP_AAS]

	# Simulate minimal media conditions throughout
	if index == N_FIXED_VARIANTS - 1:
		name = 'control_minimal'
		desc = 'Minimal media with no amino acids'
	else:
		sim_data.condition = condition_label

		# Simulate rich media conditions throughout
		if index == 0:
			sim_data.external_state.current_timeline_id = condition_label
			sim_data.external_state.saved_timelines[condition_label] = [(0, base_media_id)]

			name = 'control_rich'
			desc = 'Minimal media with amino acids control'
		# Simulate a shift from rich to minimal media
		elif index == N_FIXED_VARIANTS + len(aa_ids):
			sim_data.external_state.current_timeline_id = condition_label
			sim_data.external_state.saved_timelines[condition_label] = [
				(0, base_media_id), (SHIFT_TIME, minimal_media_id)
			]

			name = 'control_minimal_shift'
			desc = 'Remove all amino acids from rich media'
		# Simulate shifts from rich to fewer amino acids in the media
		else:
			new_media = sim_data.external_state.saved_media[minimal_media_id].copy()
			if index == 1:
				name = 'shift_to_12_aa'
				desc = 'Remove 8 AAs from rich media.'
				new_media_id = 'minimal_plus_12_amino_acids'
				for aa in AA_12:
					new_media[aa] = np.inf
			elif index == 2:
				name = 'shift_to_6_aa'
				desc = 'Remove 14 AAs from rich media.'
				new_media_id = 'minimal_plus_6_amino_acids'
				for aa in AA_6:
					new_media[aa] = np.inf
			else:
				aa = aa_ids[index - N_FIXED_VARIANTS]
				new_media_id = f'minimal_plus_{aa}'
				new_media[aa] = np.inf
				name = f'shift_to_only_{aa}'
				desc = f'Remove all but {aa} from rich media'

			sim_data.external_state.saved_media[new_media_id] = new_media

			# Create timeline to shift media at 10 minutes
			sim_data.external_state.current_timeline_id = new_media_id
			sim_data.external_state.saved_timelines[new_media_id] = [
				(0, base_media_id), (SHIFT_TIME, new_media_id)
			]

			# Add new media to lookup tables
			attrs = [
				sim_data.expectedDryMassIncreaseDict,
				sim_data.nutrient_to_doubling_time,
				sim_data.translation_supply_rate,
				sim_data.process.translation.ribosomeElongationRateDict,
				sim_data.process.translation.ribosomeFractionActiveDict,
				sim_data.process.transcription.rnaPolymeraseElongationRateDict,
				sim_data.process.transcription.rnaSynthProbFraction,
				sim_data.process.transcription.rnaSynthProbRProtein,
				sim_data.process.transcription.rnaSynthProbRnaPolymerase,
				sim_data.process.transcription.rnapFractionActiveDict,
				]
			for att in attrs:
				att[new_media_id] = att[base_media_id]

	return dict(shortName=name, desc=desc), sim_data
