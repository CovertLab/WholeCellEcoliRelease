
CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

VARIANT_OUTPUT = dict(
	shortName = "use all",
	desc = "Simulation with all kinetic constraints enabled."
	)

def kineticsConstraintsIndices(sim_data):
	return 2

def kineticsConstraints(sim_data, index):

	if index % 2 == 0:
		sim_data.process.metabolism.useAllConstraints = False
		return CONTROL_OUTPUT, sim_data
	else:
		sim_data.process.metabolism.useAllConstraints = True
		return VARIANT_OUTPUT, sim_data
