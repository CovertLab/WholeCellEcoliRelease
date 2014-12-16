
CONTROL_OUTPUT = dict(
	shortName = "wildtype",
	desc = "Wildtype simulation"
	)

def wildtypeTotalIndices(kb):
	return 0


def wildtype(kb, index):

	return dict(
		shortName = "{}_KO".format(geneID),
		desc = "Complete knockout of {}.".format(geneID)
		)
