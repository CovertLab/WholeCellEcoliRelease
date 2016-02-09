from wholecell.utils import units
from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli
from reconstruction.ecoli.fit_sim_data_1 import fitSimData_1

DOUBLING_TIMES = units.min * [100., 60., 40., 30., 24.]

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def growthRateTotalIndices(sim_data):
	return len(DOUBLING_TIMES) + 1


def growthRate(sim_data, index):
	# Refit sim_data for each growth rate

	raw_data = KnowledgeBaseEcoli()
	sim_data = fitSimData_1(raw_data, doubling_time = DOUBLING_TIMES[index])

	return dict(
		shortName = "{} min".format(DOUBLING_TIMES[index]),
		desc = "Doubling time {} min.".format(DOUBLING_TIMES[index])
		)
