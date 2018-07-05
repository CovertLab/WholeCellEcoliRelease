import collections
import numpy as np

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def tfActivityTotalIndices(sim_data):
	nNutrientTimeSeries = len(sim_data.tfToActiveInactiveConds)
	return (2 * nNutrientTimeSeries + 1)


def tfActivity(sim_data, index):

	nTfActivityTimeSeries = tfActivityTotalIndices(sim_data)

	if index % nTfActivityTimeSeries == 0:
		return CONTROL_OUTPUT, sim_data

	tfList = ["basal (no TF)"] + sorted(sim_data.tfToActiveInactiveConds)
	tf = tfList[(index + 1) // 2]
	tfStatus = None
	if index % 2 == 1:
		tfStatus = "active"
	else:
		tfStatus = "inactive"

	sim_data.condition = tf + "__" + tfStatus

	sim_data.external_state.environment.nutrients_time_series_label = tf + "__" + tfStatus
	sim_data.external_state.environment.nutrients_time_series[sim_data.external_state.environment.nutrients_time_series_label] = []
	sim_data.external_state.environment.nutrients_time_series[sim_data.external_state.environment.nutrients_time_series_label].append((
		0.0,
		sim_data.tfToActiveInactiveConds[tf][tfStatus + " nutrients"]
		))

	sim_data.genetic_perturbations = {}
	for rnaId in sim_data.conditions[sim_data.condition]["perturbations"]:
		rnaIdx = np.where(sim_data.process.transcription.rnaData["id"] == rnaId)[0]
		sim_data.genetic_perturbations[rnaId] = sim_data.process.transcription.rnaSynthProb[sim_data.condition][rnaIdx]

	return dict(
		shortName = "{}_phenotype".format(tf + "__" + tfStatus),
		desc = "Simulation of phenotype {}.".format(tf + "__" + tfStatus)
		), sim_data

