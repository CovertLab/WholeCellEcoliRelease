
def conditionIndices(sim_data):
	nConditions = len(sim_data.conditionActiveTfs)
	return nConditions

def condition(sim_data, index):

	nConditions = conditionIndices(sim_data)

	conditionLabels = sorted(sim_data.conditionActiveTfs)
	conditionLabel = conditionLabels[index]
	sim_data.condition = conditionLabel
	# TODO: add new column to condition defs to replace this?
	if sim_data.conditions[conditionLabel]["nutrients"] == "minimal_plus_amino_acids":
		sim_data.nutrientsTimeSeriesLabel = "000003_aa"

	return dict(
		shortName = "{}_env".format(conditionLabel),
		desc = "Simulation of condition {}.".format(conditionLabel)
		), sim_data
