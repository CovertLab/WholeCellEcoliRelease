
CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

from wholecell.utils import units

def EndoKcatFullRNATotalIndices(kb):
	n_variants = 3
	return n_variants + 1

def EndoKcatFullRNA(kb, index):

	# This code is not really needed but I put it here for reference
	# KcatEndoRNaseFullRNA = kb.KcatEndoRNaseFullRNA * kb.timeStep
	# KcatEndoRNaseFullRNA.checkNoUnit()
	# KcatEndoRNaseFullRNA = KcatEndoRNaseFullRNA.asNumber()

	if index == 0:
		return CONTROL_OUTPUT

	KcatEndoRNaseFullRNA = index * 0.001

	kb.KcatEndoRNaseFullRNA = KcatEndoRNaseFullRNA

	# import ipdb; ipdb.set_trace()

	# kb._parameterData['KcatEndoRNaseFullRNA'] = KcatEndoRNaseFullRNA

	kb.KcatEndoRNaseFullRNA = KcatEndoRNaseFullRNA * 1 / units.s

	return dict(
		shortName = "{} s^-1".format(KcatEndoRNaseFullRNA),
		desc = "KcatEndoRNaseFullRNA = {} s^-1.".format(KcatEndoRNaseFullRNA)
		)
