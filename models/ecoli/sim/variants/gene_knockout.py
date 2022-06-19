"""
Knockout expression of a gene

Modifies:
	sim_data.process.transcription.rna_synth_prob
	sim_data.process.transcription.rna_expression
	sim_data.process.transcription.exp_free
	sim_data.process.transcription.exp_ppgpp
	sim_data.process.transcription.attenuation_basal_prob_adjustments
	sim_data.process.transcription_regulation.basal_prob
	sim_data.process.transcription_regulation.delta_prob

Expected variant indices (depends on length of sim_data.process.transcription.rna_data):
	0: control
	1-4692: gene index to knockout
"""

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)


def gene_knockout(sim_data, index):
	rna_data = sim_data.process.transcription.rna_data

	nGenes = len(rna_data)
	nConditions = nGenes + 1

	if index % nConditions == 0:
		return CONTROL_OUTPUT, sim_data

	geneIndex = (index - 1) % nConditions
	factor = 0  # Knockout expression
	sim_data.adjust_final_expression([geneIndex], [factor])
	geneID = rna_data["id"][geneIndex]

	return dict(
		shortName = "{}_KO".format(geneID),
		desc = "Complete knockout of {}.".format(geneID)
		), sim_data
