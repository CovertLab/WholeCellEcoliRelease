"""
Variant to compare the effect of various maximum time steps.

Modifies:
	sim_data.process.replication.max_time_step
	sim_data.process.transcription.max_time_step
	sim_data.process.translation.max_time_step

Expected variant indices (dependent on TIME_STEP_FACTOR):
	0: control
	1-9: reduced max time steps (9 being lowest)
"""

TIME_STEP_FACTOR = [1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]


def time_step(sim_data, index):
	time_step_factor = TIME_STEP_FACTOR[index]
	sim_data.process.replication.max_time_step *= time_step_factor
	sim_data.process.transcription.max_time_step *= time_step_factor
	sim_data.process.translation.max_time_step *= time_step_factor

	max_time_step = min(
		sim_data.process.replication.max_time_step,
		sim_data.process.transcription.max_time_step,
		sim_data.process.translation.max_time_step,
		)

	return dict(
		shortName="{}x dt".format(time_step_factor),
		desc="Max time step adjusted by {} to {}.".format(time_step_factor, max_time_step)
		), sim_data
