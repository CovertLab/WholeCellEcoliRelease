"""
SimulationData for transcription regulation

@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import division

class TranscriptionRegulation(object):
	"""
	SimulationData for transcription regulation
	"""
	def __init__(self, raw_data, sim_data):
		# Build lookups
		self._buildLookups(raw_data, sim_data)

		# Store list of transcription factor IDs
		self.tf_ids = sim_data.tfToActiveInactiveConds.keys()

		# Build dictionary mapping transcription factors to their Kds
		self.tfKd = {}

		mRNASet = set(
			[x["id"].encode("utf-8")
			for x in raw_data.rnas
			if x["type"] != "rRNA" and x["type"] != "tRNA"]
			)

		for D in raw_data.foldChanges:
			self.tfKd[self.abbrToActiveId[D["TF"].encode("utf-8")][0]] = D["kd"]

		# Build dictionary mapping RNA targets to its regulators
		self.targetTf = {}

		for tf in sim_data.tfToFC:
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
		self.activeToBound = dict(
			[(x["active TF"].encode("utf-8"), x["metabolite bound form"].encode("utf-8"))
			for x in raw_data.tfOneComponentBound]
			)
		self.tfToTfType = dict(
			[(x["active TF"].encode("utf-8"), x["TF type"].encode("utf-8"))
			for x in raw_data.condition.tf_condition]
			)

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

	def _buildLookups(self, raw_data, sim_data):
		"""
		Builds dictionaries for mapping transcription factor abbreviations to
		their RNA IDs, and to their active form.
		"""
		geneIdToRnaId = dict(
			[(x["geneId"].encode("utf-8"), x["id"].encode("utf-8"))
			for x in raw_data.rnas]
			)

		self.abbrToRnaId = {}
		for lookupInfo in raw_data.tfIds:
			if len(lookupInfo["geneId"]) == 0:
				continue
			self.abbrToRnaId[lookupInfo["TF"].encode("utf-8")] = geneIdToRnaId[lookupInfo["geneId"].encode("utf-8")]

		self.abbrToActiveId = dict(
			[(x["TF"].encode("utf-8"), x["activeId"].encode("utf-8").split(", "))
			for x in raw_data.tfIds
			if len(x["activeId"]) > 0]
			)