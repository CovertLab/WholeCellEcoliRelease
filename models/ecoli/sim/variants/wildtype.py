
CONTROL_OUTPUT = dict(
	shortName = "wildtype",
	desc = "Wildtype simulation"
	)

def wildtypeTotalIndices(sim_data):
	return 1


def wildtype(sim_data, index):

	return CONTROL_OUTPUT, sim_data
