"""
Timelines variant for simulations with changing media at certain times

Modifies:
	sim_data.condition
	sim_data.external_state.current_timeline_id

Expected variant indices (dependent on sorted order from reconstruction/ecoli/flat/condition/timelines_def.tsv):
	0-26
"""

from __future__ import absolute_import, division, print_function


def timelines(sim_data, index):
	saved_timelines = sim_data.external_state.saved_timelines
	timeline_ids = sorted(saved_timelines)
	current_timeline_id = timeline_ids[index]
	sim_data.external_state.current_timeline_id = current_timeline_id

	# Get possible condition from starting nutrients for proper initialization
	# Not guaranteed to map to any condition or could map to multiple conditions
	nutrients = saved_timelines[current_timeline_id][0][1]
	conditions = [cond for cond in sim_data.conditionActiveTfs
		if sim_data.conditions[cond]['nutrients'] == nutrients]
	if len(conditions) == 1:
		sim_data.condition = conditions[0]
	else:
		print('Warning: could not find mapping from nutrients ({}) to condition.'
			' Using default condition ({}).'.format(nutrients, sim_data.condition))

	return dict(
		shortName = "{}_env".format(current_timeline_id),
		desc = "Simulation of environment {}.".format(current_timeline_id)
		), sim_data
