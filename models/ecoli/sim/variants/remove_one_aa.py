"""
Remove one amino acid from the minimal media plus amino acids condition.

Modifies:
	sim_data.condition
	sim_data.external_state.current_timeline_id
	sim_data.external_state.saved_timelines
	sim_data.external_state.saved_media

Expected variant indices (dependent on order of sim_data.moleculeGroups.aaIDs):
	0-20: adding one amino acid to media
	19: control (L-selenocysteine must be in media)

TODO:
	Create new media ID for new mixtures?
"""

def remove_one_aa(sim_data, index):
	# Set condition to be minimal with amino acids
	condition_label = 'with_aa'
	media_label = 'minimal_plus_amino_acids'
	sim_data.condition = condition_label
	sim_data.external_state.current_timeline_id = condition_label
	sim_data.external_state.saved_timelines[condition_label] = [(0, media_label)]

	aa = sim_data.molecule_groups.amino_acids[index][:-3]
	# Use index as a control because Sel needs to be in the media
	if aa == 'L-SELENOCYSTEINE':
		name = 'control'
		desc = 'Minimal media with amino acids control'
	# Remove one amino acid from the media
	else:
		sim_data.external_state.saved_media[media_label][aa] = 0
		name = '{}_removed'.format(aa)
		desc = 'Remove {} from rich media.'.format(aa)

	return dict(shortName=name, desc=desc), sim_data
