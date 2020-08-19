"""
SimulationData for transcription regulation

@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import absolute_import, division, print_function


class TranscriptionRegulation(object):
	"""
	SimulationData for transcription regulation
	"""
	def __init__(self, raw_data, sim_data):
		# Build lookups
		self._buildLookups(raw_data)

		# Store list of transcription factor IDs
		self.tf_ids = list(sorted(sim_data.tfToActiveInactiveConds.keys()))

		# Build dictionary mapping transcription factors to their Kds
		self.tfKd = {}

		mRNASet = {
			x["id"]
			for x in raw_data.rnas
			if x["type"] not in ("rRNA", "tRNA")}

		for D in raw_data.foldChanges:
			self.tfKd[self.abbrToActiveId[D["TF"]][0]] = D["kd"]

		# Build dictionary mapping RNA targets to its regulators
		self.targetTf = {}

		for tf in sorted(sim_data.tfToFC):
			targets = sim_data.tfToFC[tf]
			targetsToRemove = []

			for target in targets:
				if target not in mRNASet:
					targetsToRemove.append(target)
					continue

				if target not in self.targetTf:
					self.targetTf[target] = []

				self.targetTf[target].append(tf)

			for targetToRemove in targetsToRemove:
				sim_data.tfToFC[tf].pop(targetToRemove)

		# Build dictionaries mapping transcription factors to their bound form,
		# and to their regulating type
		self.activeToBound = {
			x["active TF"]: x["metabolite bound form"]
			for x in raw_data.tfOneComponentBound}
		self.tfToTfType = {
			x["active TF"]: x["TF type"]
			for x in raw_data.condition.tf_condition}

	def pPromoterBoundTF(self, tfActive, tfInactive):
		"""
		Computes probability of a transcription factor binding promoter.
		"""
		return float(tfActive) / (float(tfActive) + float(tfInactive))

	def pPromoterBoundSKd(self, signal, Kd, power):
		"""
		Computes probability of a one-component transcription factor binding
		promoter.
		"""
		return float(signal)**power / (float(signal)**power + float(Kd))

	def _buildLookups(self, raw_data):
		"""
		Builds dictionaries for mapping transcription factor abbreviations to
		their RNA IDs, and to their active form.
		"""
		geneIdToRnaId = {x["id"]: x["rnaId"] for x in raw_data.genes}

		self.abbrToRnaId = {}
		for lookupInfo in raw_data.tfIds:
			if len(lookupInfo["geneId"]) == 0:
				continue
			self.abbrToRnaId[lookupInfo["TF"]] = geneIdToRnaId[lookupInfo["geneId"]]

		self.abbrToActiveId = {
			x["TF"]: x["activeId"].split(", ")
			for x in raw_data.tfIds
			if len(x["activeId"]) > 0}
