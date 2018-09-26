"""
Variant to compare impact of adjusting the kinetics weighting factor in metabolism.
Weighting varies between 0 (no kinetics) and 1 (only kinetics).

Modifies:
	sim_data.process.metabolism.kinetic_objective_weight

Expected variant indices (dependent on KINETIC_OBJECTIVE_WEIGHT):
	3: control
	0 (no kinetics) - 9 (only kinetics)
"""

KINETIC_OBJECTIVE_WEIGHT = [0, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 0.01, 0.1, 1]


def metabolism_kinetic_objective_weight(sim_data, index):
	weight = KINETIC_OBJECTIVE_WEIGHT[index]
	sim_data.process.metabolism.kinetic_objective_weight = weight

	return dict(
		shortName="lambda={:.0E}".format(weight),
		desc="Simulation with kinetics objective weight of {}.".format(weight)
		), sim_data
