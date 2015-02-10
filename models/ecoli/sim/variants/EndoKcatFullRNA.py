
CONTROL_OUTPUT = dict(
	shortName = "EndoKcatFullRNA",
	desc = "Sensitivity analysis Kcat"
	)

def EndoKcatFullRNATotalIndices(kb):
	return 1

def EndoKcatFullRNA(kb, index):

	KcatEndoRNaseFullRNA = index * 0.001


	kb.KcatEndoRNaseFullRNA = KcatEndoRNaseFullRNA

	return CONTROL_OUTPUT