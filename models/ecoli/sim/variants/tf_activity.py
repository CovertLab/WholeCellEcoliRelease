"""
Variant for analyzing transcription factor activity in active and inactive
conditions. Will set the condition for each TF and the appropriate nutrients
and/or genetic perturbations for the TF to be active or inactive.

Modifies:
	sim_data.condition
	sim_data.external_state.current_timeline_id
	sim_data.external_state.saved_timelines
	sim_data.genetic_perturbations

Expected variant indices (dependent on length of sim_data.tfToActiveInactiveConds):
	0: control
	1+: modify condition for one transcription factor
		(odd values will be active TF, even values will be inactive TF)
"""

from __future__ import absolute_import, division, print_function

import numpy as np

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)


def tf_activity(sim_data, index):
	nNutrientTimeSeries = len(sim_data.tfToActiveInactiveConds)
	nTfActivityTimeSeries = (2 * nNutrientTimeSeries + 1)

	if index % nTfActivityTimeSeries == 0:
		return CONTROL_OUTPUT, sim_data

	tfList = ["basal (no TF)"] + sorted(sim_data.tfToActiveInactiveConds)
	tf = tfList[(index + 1) // 2]
	if index % 2 == 1:
		tfStatus = "active"
	else:
		tfStatus = "inactive"

	sim_data.condition = tf + "__" + tfStatus

	sim_data.external_state.current_timeline_id = tf + "__" + tfStatus
	sim_data.external_state.saved_timelines[sim_data.external_state.environment.current_timeline_id] = []
	sim_data.external_state.saved_timelines[sim_data.external_state.environment.current_timeline_id].append((
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

