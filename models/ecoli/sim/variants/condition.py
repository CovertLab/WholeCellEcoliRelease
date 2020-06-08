"""
Condition variant for simulations in different environmental conditions

Modifies:
	sim_data.condition
	sim_data.external_state.current_timeline_id

Expected variant indices (dependent on sim_data.ordered_conditions and should
be the same order as rows in condition_defs.tsv):
	0: control
	1: anaerobic
	2: with amino acids
	3: acetate
	4: succinate
"""
from __future__ import absolute_import, division, print_function


def condition(sim_data, index):
	condition_labels = sim_data.ordered_conditions
	condition_label = condition_labels[index]
	sim_data.condition = condition_label
	sim_data.external_state.current_timeline_id = condition_label
	sim_data.external_state.saved_timelines[condition_label] = [
		(0, sim_data.conditions[condition_label]["nutrients"])
	]

	return dict(
		shortName = "{}_env".format(condition_label),
		desc = "Simulation of condition {}.".format(condition_label)
		), sim_data
