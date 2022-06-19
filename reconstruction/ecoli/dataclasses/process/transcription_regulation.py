"""
SimulationData for transcription regulation

"""

from typing import Union

import numpy as np
from scipy import sparse


class TranscriptionRegulation(object):
	"""
	SimulationData for transcription regulation
	"""
	def __init__(self, raw_data, sim_data):
		# Build lookups
		self._build_lookups(raw_data)

		# Store list of transcription factor IDs
		self.tf_ids = list(sorted(sim_data.tf_to_active_inactive_conditions.keys()))

		# Build dictionary mapping RNA targets to its regulators
		self.target_tf = {}

		for tf in sorted(sim_data.tf_to_fold_change):
			targets = sim_data.tf_to_fold_change[tf]
			targetsToRemove = []

			for target in targets:
				if target not in self.target_tf:
					self.target_tf[target] = []

				self.target_tf[target].append(tf)

			for targetToRemove in targetsToRemove:
				sim_data.tf_to_fold_change[tf].pop(targetToRemove)

		# Build dictionaries mapping transcription factors to their bound form,
		# and to their regulating type
		self.active_to_bound = {
			x["active TF"]: x["metabolite bound form"]
			for x in raw_data.tf_one_component_bound}
		self.tf_to_tf_type = {
			x["active TF"]: x["TF type"]
			for x in raw_data.condition.tf_condition}
		self.tf_to_gene_id = {
			x["active TF"]: x["TF"]
			for x in raw_data.condition.tf_condition}

		# Values set after promoter fitting in parca with calculateRnapRecruitment()
		self.basal_prob = None
		self.delta_prob = None

	def p_promoter_bound_tf(self, tfActive, tfInactive):
		"""
		Computes probability of a transcription factor binding promoter.
		"""
		return float(tfActive) / (float(tfActive) + float(tfInactive))

	def p_promoter_bound_SKd(self, signal, Kd, power):
		"""
		Computes probability of a one-component transcription factor binding
		promoter.
		"""
		return float(signal)**power / (float(signal)**power + float(Kd)**power)

	def get_delta_prob_matrix(self, dense=False, ppgpp=False) -> Union[sparse.csr_matrix, np.ndarray]:
		"""
		Returns the delta probability matrix mapping the promoter binding effect
		of each TF to each gene.

		Args:
			dense: If True, returns a dense matrix, otherwise csr sparse
			ppgpp: If True, normalizes delta probabilities to be on the same
				scale as ppGpp normalized probabilities since delta_prob is
				calculated based on basal_prob which is not normalized to 1

		Returns:
			delta_prob: matrix of probabilities changes expected with a TF
				binding to a promoter for each gene (n genes, m TFs)
		"""

		ppgpp_scaling = self.basal_prob[self.delta_prob['deltaI']]
		ppgpp_scaling[ppgpp_scaling == 0] = 1
		scaling_factor = ppgpp_scaling if ppgpp else 1.
		delta_prob = sparse.csr_matrix(
			(self.delta_prob['deltaV'] / scaling_factor,
			(self.delta_prob['deltaI'], self.delta_prob['deltaJ'])),
			shape=self.delta_prob['shape'])

		if dense:
			delta_prob = delta_prob.toarray()

		return delta_prob

	def _build_lookups(self, raw_data):
		"""
		Builds dictionaries for mapping transcription factor abbreviations to
		their RNA IDs, and to their active form.
		"""
		gene_id_to_cistron_id = {x['id']: x['rna_ids'][0] for x in raw_data.genes}

		self.abbr_to_rna_id = {}
		for lookupInfo in raw_data.transcription_factors:
			if len(lookupInfo["geneId"]) == 0 or lookupInfo["geneId"] not in gene_id_to_cistron_id:
				continue
			self.abbr_to_rna_id[lookupInfo["TF"]] = gene_id_to_cistron_id[lookupInfo["geneId"]]

		self.abbr_to_active_id = {
			x["TF"]: x["activeId"].split(", ")
			for x in raw_data.transcription_factors
			if len(x["activeId"]) > 0}
