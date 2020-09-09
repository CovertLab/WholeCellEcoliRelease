"""
SimulationData for transcription regulation

"""

from __future__ import absolute_import, division, print_function


class TranscriptionRegulation(object):
	"""
	SimulationData for transcription regulation
	"""
	def __init__(self, raw_data, sim_data):
		# Build lookups
		self._build_lookups(raw_data)

		# Store list of transcription factor IDs
		self.tf_ids = list(sorted(sim_data.tf_to_active_inactive_conditions.keys()))

		# Build dictionary mapping transcription factors to their Kds
		self.tf_Kd = {}

		mRNASet = {
			x["id"]
			for x in raw_data.rnas
			if x["type"] not in ("rRNA", "tRNA")}

		for D in raw_data.fold_changes:
			self.tf_Kd[self.abbr_to_active_id[D["TF"]][0]] = D["kd"]

		# Build dictionary mapping RNA targets to its regulators
		self.target_tf = {}

		for tf in sorted(sim_data.tf_to_fold_change):
			targets = sim_data.tf_to_fold_change[tf]
			targetsToRemove = []

			for target in targets:
				if target not in mRNASet:
					targetsToRemove.append(target)
					continue

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
		return float(signal)**power / (float(signal)**power + float(Kd))

	def _build_lookups(self, raw_data):
		"""
		Builds dictionaries for mapping transcription factor abbreviations to
		their RNA IDs, and to their active form.
		"""
		geneIdToRnaId = {x["id"]: x['rna_id'] for x in raw_data.genes}

		self.abbr_to_rna_id = {}
		for lookupInfo in raw_data.transcription_factors:
			if len(lookupInfo["geneId"]) == 0:
				continue
			self.abbr_to_rna_id[lookupInfo["TF"]] = geneIdToRnaId[lookupInfo["geneId"]]

		self.abbr_to_active_id = {
			x["TF"]: x["activeId"].split(", ")
			for x in raw_data.transcription_factors
			if len(x["activeId"]) > 0}
