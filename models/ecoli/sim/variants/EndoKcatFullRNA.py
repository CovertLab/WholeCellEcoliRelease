
CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

from wholecell.utils import units

def EndoKcatFullRNATotalIndices(sim_data):
	n_variants = 3
	return n_variants + 1

def EndoKcatFullRNA(sim_data, index):

	# This code is not really needed but I put it here for reference
	# KcatEndoRNaseFullRNA = sim_data.KcatEndoRNaseFullRNA * sim_data.timeStep
	# KcatEndoRNaseFullRNA.checkNoUnit()
	# KcatEndoRNaseFullRNA = KcatEndoRNaseFullRNA.asNumber()

	if index == 0:
		return CONTROL_OUTPUT

	KcatEndoRNaseFullRNA = index * 0.001

	sim_data.KcatEndoRNaseFullRNA = KcatEndoRNaseFullRNA

	# import ipdb; ipdb.set_trace()

	# sim_data._parameterData['KcatEndoRNaseFullRNA'] = KcatEndoRNaseFullRNA

	sim_data.KcatEndoRNaseFullRNA = KcatEndoRNaseFullRNA * 1 / units.s

	return dict(
		shortName = "{} s^-1".format(KcatEndoRNaseFullRNA),
		desc = "KcatEndoRNaseFullRNA = {} s^-1.".format(KcatEndoRNaseFullRNA)
		)
