CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)


def nutrient_time_series(sim_data, index):
	n_nutrients_time_series = len(sim_data.external_state.environment.nutrients_time_series)

	if index % n_nutrients_time_series == 0:
		return CONTROL_OUTPUT, sim_data

	nutrients_time_series_labels = sorted(sim_data.external_state.environment.nutrients_time_series)
	nutrients_time_series_label = nutrients_time_series_labels[index]
	sim_data.external_state.environment.nutrients_time_series_label = nutrients_time_series_label

	return dict(
		shortName = "{}_env".format(nutrients_time_series_label),
		desc = "Simulation of environment {}.".format(nutrients_time_series_label)
		), sim_data
