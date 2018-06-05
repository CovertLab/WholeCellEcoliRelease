
def conditionIndices(sim_data):
	nConditions = len(sim_data.conditionActiveTfs)
	return nConditions

def condition(sim_data, index):

	nConditions = conditionIndices(sim_data)

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
