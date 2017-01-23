CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def nutrientTimeSeriesDownshiftTotalIndices(sim_data):
	nNutrientTimeSeries = len(sim_data.nutrientsTimeSeries)
	return nNutrientTimeSeries + 1


def nutrientTimeSeriesDownshift(sim_data, index):

	nNutrientTimeSeries = nutrientTimeSeriesDownshiftTotalIndices(sim_data)

	if index % nNutrientTimeSeries == 0:
		return CONTROL_OUTPUT, sim_data

	# Start in minimal_plus_amino acids, index = 2
	sim_data.condition = sorted(sim_data.conditionActiveTfs)[2]

	nutrientTimeSeriesLabels = sorted(sim_data.nutrientsTimeSeries)
	nutrientTimeSeriesLabel = nutrientTimeSeriesLabels[index]
	sim_data.nutrientsTimeSeriesLabel = nutrientTimeSeriesLabel

	return dict(
		shortName = "{}_env".format(nutrientTimeSeriesLabel),
		desc = "Simulation of environment {}.".format(nutrientTimeSeriesLabel)
		), sim_data
