
CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def environmentTotalIndices(sim_data):
	nEnvironments = len(sim_data.envDict)
	return nEnvironments


def environment(sim_data, index):
	# Knocks-out genes in order

	nEnvironments = environmentTotalIndices(sim_data)

	if index % nEnvironments == 0:
		return CONTROL_OUTPUT, sim_data

	environmentNames = sorted(sim_data.envDict)
	envName = environmentNames[index]
	sim_data.environment = envName

	return dict(
		shortName = "{}_env".format(envName),
		desc = "Simulation of environment {}.".format(envName)
		), sim_data
