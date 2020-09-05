"""Utilities for metrics pickle file

Computing metrics requires some data from sim_data, but we don't want to
load all of sim_data when computing metrics. Instead, we save off a
subset of sim_data into a metrics pickle file. These functions are used
to generate that file.
"""

from __future__ import absolute_import, division, print_function


def get_metrics_data_dict(sim_data):
	"""Create a dictionary to save as the metrics data pickle.

	Arguments:
		sim_data: The sim_data object the dictionary values will be
			drawn from.

	Returns:
		Dictionary suitable for serializing as metrics data pickle.
	"""
	is_mRNA = sim_data.process.transcription.rna_data['is_mRNA']
	metrics_data = {
		"translation_monomer_ids":
			sim_data.process.translation.monomer_data["id"],
		"rna_ids":
			sim_data.process.transcription.rna_data["id"][is_mRNA],
		"expected_mRNA_counts":
			sim_data.process.transcription.rna_expression[sim_data.condition][is_mRNA],
	}
	return metrics_data
