"""
Remove one amino acid to the minimal media plus amino acids condition after 10 minutes.

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

Expected variant indices (dependent on order of sim_data.moleculeGroups.aaIDs):
	0-20: adding one amino acid to media
	19: control (adding L-selenocysteine which is already required in media)
"""


SHIFT_TIME = 600  # 10 minutes


def remove_one_aa_shift(sim_data, index):
	# Remove one amino acid from the media
	condition_label = 'with_aa'
	sim_data.condition = condition_label
	base_media_id = 'minimal_plus_amino_acids'

	aa = sim_data.molecule_groups.amino_acids[index][:-3]
	# Use index as a control because Sel needs to be in the media
	if aa == 'L-SELENOCYSTEINE':
		sim_data.external_state.current_timeline_id = condition_label
		sim_data.external_state.saved_timelines[condition_label] = [(0, base_media_id)]

		name = 'control'
		desc = 'Minimal media with amino acids control'
	# Remove one amino acid from the media
	else:
		new_media_id = f'{base_media_id}_plus_{aa}'
		new_media = sim_data.external_state.saved_media[base_media_id].copy()
		new_media[aa] = 0
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

		name = 'shift_{}_removed'.format(aa)
		desc = 'Remove {} from rich media.'.format(aa)

	return dict(shortName=name, desc=desc), sim_data
