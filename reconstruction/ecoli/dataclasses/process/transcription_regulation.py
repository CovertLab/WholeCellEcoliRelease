
from __future__ import division

import numpy as np
import os
import cPickle
import wholecell
from wholecell.utils import units
import collections

class TranscriptionRegulation(object):
	def __init__(self, raw_data, sim_data):
		self._buildLookups(raw_data, sim_data)

		# self.tfs = sorted(set([x["TF"].encode("utf-8") for x in raw_data.foldChanges]))
		# self.tfTarget = collections.defaultdict(list)
		# self.targetTF = collections.defaultdict(list)
		# self.tfIdTargetId = collections.defaultdict(list)
		# self.targetIdTFId = collections.defaultdict(list)
		self.tfKd = {}
		for D in raw_data.foldChanges:
			# self.tfTarget[ D["TF"].encode("utf-8") ].append( D["Target"].encode("utf-8") )
			# self.targetTF[ D["Target"].encode("utf-8") ].append( D["TF"].encode("utf-8") )

			# # Ignore things like two-component systems
			# if D["TF"].encode("utf-8") not in self.abbrToActiveId:
			# 	continue

			# for activeId in self.abbrToActiveId[ D["TF"].encode("utf-8") ]:
			# 	self.tfIdTargetId[activeId].append(self.abbrToRnaId[ D["Target"].encode("utf-8") ])
			# 	self.targetIdTFId[self.abbrToRnaId[ D["Target"].encode("utf-8") ]].append(activeId)

			self.tfKd[self.abbrToActiveId[D["TF"].encode("utf-8")][0]] = D["kd"]

		self.targetTf = {}
		for tf in sim_data.tfToFC:
			targets = sim_data.tfToFC[tf]
			for target in targets:
				if target not in self.targetTf:
					self.targetTf[target] = []
				self.targetTf[target].append(tf)

		self.tfNTargets = dict([(key, len(val)) for key,val in sim_data.tfToFC.iteritems()])
		self.activeToBound = dict([(x["active TF"].encode("utf-8"), x["metabolite bound form"].encode("utf-8")) for x in raw_data.tfOneComponentBound])
		return

	def pPromoterBoundTF(self, tfActive, tfInactive):
		return float(tfActive) / (float(tfActive) + float(tfInactive))

	def pPromoterBoundSKd(self, signal, Kd, power):
		return float(signal)**power / (float(signal)**power + float(Kd))

	def pPromoterBound(self, dissocConstant, nPromoters, nTfs):
		b = float(dissocConstant + nPromoters + nTfs)
		c = float(nPromoters * nTfs)
		pos = (b + np.sqrt(b**2 - 4 * c)) / (2 * nPromoters)
		neg = (b - np.sqrt(b**2 - 4 * c)) / (2 * nPromoters)
		if (0 <= pos and pos <= 1) and (0 <= neg and neg <= 1):
			raise Exception, "Both solutions to the quadratic equation are between 0 and 1"
		if (0 <= pos and pos <= 1):
			return pos
		if (0 <= neg and neg <= 1):
			return neg
		raise Exception, "Neither solution to the quadratic equation is between 0 and 1"

	def pTfBound(self, dissocConstant, nPromoters, nTfs):
		b = float(dissocConstant + nPromoters + nTfs)
		c = float(nPromoters * nTfs)
		pos = (b + np.sqrt(b**2 - 4 * c)) / (2 * nTfs)
		neg = (b - np.sqrt(b**2 - 4 * c)) / (2 * nTfs)
		if (0 <= pos and pos <= 1) and (0 <= neg and neg <= 1):
			raise Exception, "Both solutions to the quadratic equation are between 0 and 1"
		if (0 <= pos and pos <= 1):
			return pos
		if (0 <= neg and neg <= 1):
			return neg
		raise Exception, "Neither solution to the quadratic equation is between 0 and 1"

	def _buildLookups(self, raw_data, sim_data):
		geneIdToRnaId = dict([(x["geneId"].encode("utf-8"), x["id"].encode("utf-8")) for x in raw_data.rnas])

		self.abbrToRnaId = {}
		for lookupInfo in raw_data.tfIds:
			if len(lookupInfo["geneId"]) == 0:
				continue
			self.abbrToRnaId[lookupInfo["TF"].encode("utf-8")] = geneIdToRnaId[lookupInfo["geneId"].encode("utf-8")]

		self.abbrToActiveId = dict([(x["TF"].encode("utf-8"), x["activeId"].encode("utf-8").split(", ")) for x in raw_data.tfIds if len(x["activeId"]) > 0])

		# # Get rid of two-component systems, which we will incorporate later
		# keysToDelete = []
		# excludes = ["MONOMER0-1",
		# 			"PROTEIN-NRIP",
		# 			"MONOMER0-4180",
		# 			"CPLX0-7754",
		# 			"CPLX0-7884",
		# 			"MONOMER0-4276",
		# 			"CPLX0-7748",
		# 			"MONOMER0-4288",
		# 			"CPLX0-7721",
		# 			"MONOMER0-2741",
		# 			"MONOMER0-4147",
		# 			"CPLX0-7795",
		# 			"MONOMER0-4341",
		# 			]
		# for key, val in self.abbrToActiveId.iteritems():
		# 	if val[0].startswith("PHOSPHO") or val[0] in excludes:
		# 		keysToDelete.append(key)
		# 		continue
		# 	locations = sim_data.getter.getLocation(val)
		# 	L = []
		# 	for frameId, loc in zip(val, locations):
		# 		L.append(frameId + "[" + loc[0] + "]")
		# 	self.abbrToActiveId[key] = L

		# for key in keysToDelete:
		# 	del self.abbrToActiveId[key]