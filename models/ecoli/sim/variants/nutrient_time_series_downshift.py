CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def nutrientTimeSeriesDownshiftTotalIndices(sim_data):
	nNutrientTimeSeries = len(sim_data.external_state.environment.nutrients_time_series)
	return nNutrientTimeSeries + 1


def nutrientTimeSeriesDownshift(sim_data, index):

	nNutrientTimeSeries = nutrientTimeSeriesDownshiftTotalIndices(sim_data)

	if index % nNutrientTimeSeries == 0:
		return CONTROL_OUTPUT, sim_data

	# Start in minimal_plus_amino acids, index = 2
	sim_data.condition = sorted(sim_data.conditionActiveTfs)[2]

	nutrients_time_series_labels = sorted(sim_data.external_state.environment.nutrients_time_series)
	nutrients_time_series_label = nutrients_time_series_labels[index]
	sim_data.external_state.environment.nutrients_time_series_label = nutrients_time_series_label

	return dict(
		shortName = "{}_env".format(nutrients_time_series_label),
		desc = "Simulation of environment {}.".format(nutrients_time_series_label)
		), sim_data
