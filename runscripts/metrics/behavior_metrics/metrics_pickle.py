def gen_metrics_data_dict(sim_data):
	is_mRNA = sim_data.process.transcription.rnaData["isMRna"]
	metrics_data = {
		"translation_monomer_ids":
			sim_data.process.translation.monomerData["id"],
		"rna_ids":
			sim_data.process.transcription.rnaData["id"][is_mRNA],
		"expected_mRNA_counts":
			sim_data.process.transcription.rnaExpression[sim_data.condition][is_mRNA],
	}
	return metrics_data
