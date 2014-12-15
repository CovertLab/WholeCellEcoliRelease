
CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def geneKnockout(kb, index):
	# Knocks-out genes in order

	nGenes = kb.rnaData.fullArray().size

	nConditions = nGenes + 1

	if (index + nGenes) % nConditions == 0:
		return CONTROL_OUTPUT

	geneIndex = (index - 1) % nConditions

	synthProb = kb.rnaData["synthProb"]

	synthProb[geneIndex] = 0.0

	synthProb /= synthProb.sum()

	kb.rnaData["synthProb"] = synthProb

	geneID = kb.rnaData["id"][geneIndex]

	return dict(
		shortName = "{}_KO".format(geneID),
		desc = "Complete knockout of {}.".format(geneID)
		)
