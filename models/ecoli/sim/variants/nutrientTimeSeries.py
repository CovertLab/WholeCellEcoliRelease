
CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def nutrientTimeSeriesTotalIndices(sim_data):
	nNutrientTimeSeries = len(sim_data.nutrientsTimeSeries)
	return nNutrientTimeSeries


def nutrientTimeSeries(sim_data, index):

	nNutrientTimeSeries = nutrientTimeSeriesTotalIndices(sim_data)

	if index % nNutrientTimeSeries == 0:
		return CONTROL_OUTPUT, sim_data

	nutrientTimeSeriesLabels = sorted(sim_data.nutrientsTimeSeries)
	nutrientTimeSeriesLabel = nutrientTimeSeriesLabels[index]
	sim_data.nutrientsTimeSeriesLabel = nutrientTimeSeriesLabel

	return dict(
		shortName = "{}_env".format(nutrientTimeSeriesLabel),
		desc = "Simulation of environment {}.".format(nutrientTimeSeriesLabel)
		), sim_data
