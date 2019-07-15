"""
Condition variant for simulations in different environmental conditions

Modifies:
	sim_data.condition
	sim_data.external_state.environment.nutrients_time_series_label

Expected variant indices (dependent on sorted order of sim_data.conditionActiveTfs):
	0: control
	1: anaerobic
	2: with amino acids
"""

def conditionIndices(sim_data):
	nConditions = len(sim_data.conditionActiveTfs)
	return nConditions

def condition(sim_data, index):
	condition_labels = sorted(sim_data.conditionActiveTfs)
	condition_label = condition_labels[index]
	sim_data.condition = condition_label
	# TODO: add new column to condition defs to replace this?
	if sim_data.conditions[condition_label]["nutrients"] == "minimal_plus_amino_acids":
		sim_data.external_state.environment.nutrients_time_series_label = "000003_aa"
	elif sim_data.conditions[condition_label]["nutrients"] == "minimal_minus_oxygen":
		sim_data.external_state.environment.nutrients_time_series_label = "000004_oxygen_absent"

	return dict(
		shortName = "{}_env".format(condition_label),
		desc = "Simulation of condition {}.".format(condition_label)
		), sim_data
