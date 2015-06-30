TIME_STEPS = [0.1, 0.25, 0.5, 0.75]

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def timeStepTotalIndices(kb):
	nTimeSteps = len(TIME_STEPS)
	nConditions = nTimeSteps + 1
	return nConditions


def timeStep(kb, index):
	# Vary time step used for simulation

	nConditions = timeStepTotalIndices(kb)

	if index == 0:
		return CONTROL_OUTPUT

	kb.timeStepSec = TIME_STEPS[index - 1]

	return dict(
		shortName = "{} sec".format(TIME_STEPS[index - 1]),
		desc = "Simulation uses time step of {} seconds.".format(TIME_STEPS[index - 1])
		)
