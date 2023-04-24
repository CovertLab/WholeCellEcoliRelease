"""
Amino acid cut using different trna synthetase kinetic solutions

Modifies:
	sim_data.condition
	sim_data.external_state.current_timeline_id

"""

from __future__ import absolute_import, division, print_function


index_to_sweep_level = {
		0: 2,
		1: 3,
		2: 4,
		3: 5,
	}

def trna_synthetase_kinetics_test(sim_data, index):
	# Set timeline to amino acid cut (index 25)
	saved_timelines = sim_data.external_state.saved_timelines
	timeline_ids = sorted(saved_timelines)
	current_timeline_id = timeline_ids[25]
	sim_data.external_state.current_timeline_id = current_timeline_id

	nutrients = saved_timelines[current_timeline_id][0][1]
	conditions = [cond for cond in sim_data.condition_active_tfs
		if sim_data.conditions[cond]['nutrients'] == nutrients]
	sim_data.condition = conditions[0]

	# Set trna synthetase kinetics
	sweep_level = index_to_sweep_level[index]
	kinetics = sim_data.relation.trna_charging_kinetics

	k_cats_dict = kinetics[sweep_level]['synthetase_to_k_cat']
	sim_data.relation.synthetase_to_k_cat = k_cats_dict

	K_M_As_dict = kinetics[sweep_level]['synthetase_to_K_A']
	sim_data.relation.synthetase_to_K_A = K_M_As_dict

	K_M_Ts_dict = kinetics[sweep_level]['trna_to_K_T']
	sim_data.relation.trna_to_K_T = K_M_Ts_dict

	short_name = f'Testing trna synthetase kinetics, sweep_level = {sweep_level}, aa cut'

	return dict(
		shortName=short_name,
		desc=f'tRNA Synthetase Kinetics Test Variant is set at index {index}.'
		), sim_data
