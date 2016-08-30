import collections

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

	sim_data.nutrientsTimeSeriesLabel = tf + "__" + tfStatus
	sim_data.nutrientsTimeSeries[sim_data.nutrientsTimeSeriesLabel] = collections.deque()
	sim_data.nutrientsTimeSeries[sim_data.nutrientsTimeSeriesLabel].append((
		0.0,
		sim_data.tfToActiveInactiveConds[tf][tfStatus + " nutrients"]
		))

	return dict(
		shortName = "{}_phenotype".format(tf + "__" + tfStatus),
		desc = "Simulation of phenotype {}.".format(tf + "__" + tfStatus)
		), sim_data

