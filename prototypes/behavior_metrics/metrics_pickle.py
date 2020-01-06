def gen_metrics_data_dict(sim_data):
	metrics_data = {
		"translation_monomer_ids": sim_data.process.translation.monomerData["id"]
	}
	return metrics_data
