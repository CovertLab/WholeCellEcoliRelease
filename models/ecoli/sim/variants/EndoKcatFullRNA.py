
CONTROL_OUTPUT = dict(
	shortName = "EndoKcatFullRNA",
	desc = "Sensitivity analysis Kcat"
	)

from wholecell.utils import units

def EndoKcatFullRNATotalIndices(kb):
	return 1

def EndoKcatFullRNA(kb, index):

	#KcatEndoRNaseFullRNA = index * 0.001

	#kb._parameterData["KcatEndoRNaseFullRNA"] = KcatEndoRNaseFullRNA
	
	timeStepSec = 1 # TODO

	KcatEndoRNaseFullRNA = kb.KcatEndoRNaseFullRNA.asNumber(1 / units.s) * timeStepSec
	KcatEndoRNaseFullRNA = index * 0.001
	kb.KcatEndoRNaseFullRNA = KcatEndoRNaseFullRNA

	# import ipdb; ipdb.set_trace()

	# kb._parameterData['KcatEndoRNaseFullRNA'] = KcatEndoRNaseFullRNA

	return CONTROL_OUTPUT