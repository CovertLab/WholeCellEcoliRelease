TIME_STEPS = [0.25, 0.5, 0.75]

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def timeStepTotalIndices(sim_data):
	nTimeSteps = len(TIME_STEPS)
	nConditions = nTimeSteps + 1
	return nConditions


def timeStep(sim_data, index):
	# Vary time step used for simulation

	nConditions = timeStepTotalIndices(sim_data)

	if index == 0:
		return CONTROL_OUTPUT, sim_data

	sim_data.timeStepSec = TIME_STEPS[index - 1]

	return dict(
		shortName = "{} sec".format(TIME_STEPS[index - 1]),
		desc = "Simulation uses time step of {} seconds.".format(TIME_STEPS[index - 1])
		), sim_data
