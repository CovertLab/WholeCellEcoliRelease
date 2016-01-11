
CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

from wholecell.utils import units

def EndoKcatFullRNATotalIndices(sim_data):
	n_variants = 3
	return n_variants + 1

def EndoKcatFullRNA(sim_data, index):

	if index == 0:
		return CONTROL_OUTPUT

	KcatEndoRNaseFullRNA = index * 0.001
	sim_data.KcatEndoRNaseFullRNA = KcatEndoRNaseFullRNA * 1 / units.s

	return dict(
		shortName = "{} s^-1".format(KcatEndoRNaseFullRNA),
		desc = "KcatEndoRNaseFullRNA = {} s^-1.".format(KcatEndoRNaseFullRNA)
		)
