CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)


def timelines(sim_data, index):
	'''
	Run variants of different timelines, saved in condition/media/timelines_def.tsv
	Args:
		index: the index in timelines_def.tsv for this variant's timeline

	'''

	n_saved_timelines = len(sim_data.external_state.environment.saved_timelines)

	if index % n_saved_timelines == 0:
		return CONTROL_OUTPUT, sim_data

	timeline_ids = sorted(sim_data.external_state.environment.saved_timelines)
	current_timeline_id = timeline_ids[index]
	sim_data.external_state.environment.current_timeline_id = current_timeline_id

	return dict(
		shortName = "{}_env".format(current_timeline_id),
		desc = "Simulation of environment {}.".format(current_timeline_id)
		), sim_data
