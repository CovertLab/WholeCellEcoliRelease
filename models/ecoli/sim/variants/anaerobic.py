
AEROBIC_OUTPUT = dict(
	shortName = "aerobic",
	desc = "Aerobic simulation"
	)

ANAEROBIC_OUTPUT = dict(
	shortName = "anaerobic",
	desc = "Anaerobic simulation"
	)

def anaerobicIndices(sim_data):
	return 2

def anaerobic(sim_data, index):

	if index % 2 == 0:
		return AEROBIC_OUTPUT, sim_data

	sim_data.nutrientsTimeSeriesLabel = "000004_oxygen_absent"
	sim_data.condition = "no_oxygen"
	sim_data.doubling_time = sim_data.conditionToDoublingTime[sim_data.condition]

	return ANAEROBIC_OUTPUT, sim_data
