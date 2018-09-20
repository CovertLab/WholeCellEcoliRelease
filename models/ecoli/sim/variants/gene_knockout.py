CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)


def gene_knockout(sim_data, index):
	# Knocks-out genes in order
	nGenes = sim_data.process.transcription.rnaData.fullArray().size
	nConditions = nGenes + 1

	if index % nConditions == 0:
		return CONTROL_OUTPUT, sim_data

	geneIndex = (index - 1) % nConditions

	synthProb = sim_data.process.transcription.rnaData["synthProb"]

	synthProb[geneIndex] = 0.0

	synthProb /= synthProb.sum()

	sim_data.process.transcription.rnaData["synthProb"] = synthProb

	geneID = sim_data.process.transcription.rnaData["id"][geneIndex]

	return dict(
		shortName = "{}_KO".format(geneID),
		desc = "Complete knockout of {}.".format(geneID)
		), sim_data
