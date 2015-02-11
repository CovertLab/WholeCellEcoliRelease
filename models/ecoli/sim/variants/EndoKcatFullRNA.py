
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
	
	KcatEndoRNaseFullRNA = kb.KcatEndoRNaseFullRNA.asNumber(1 / units.s) * self.timeStepSec
	KcatEndoRNaseFullRNA = index * 0.001

	return CONTROL_OUTPUT