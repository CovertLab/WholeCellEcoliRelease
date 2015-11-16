
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

		self.tfs = sorted(set([x["TF"].encode("utf-8") for x in raw_data.foldChanges]))
		self.tfTarget = collections.defaultdict(list)
		self.targetTF = collections.defaultdict(list)
		self.tfIdTargetId = collections.defaultdict(list)
		self.targetIdTFId = collections.defaultdict(list)
		self.tfKd = {}
		for D in raw_data.foldChanges:
			self.tfTarget[ D["TF"].encode("utf-8") ].append( D["Target"].encode("utf-8") )
			self.targetTF[ D["Target"].encode("utf-8") ].append( D["TF"].encode("utf-8") )

			# Ignore things like two-component systems
			if D["TF"].encode("utf-8") not in self.abbrToActiveId:
				continue

			for activeId in self.abbrToActiveId[ D["TF"].encode("utf-8") ]:
				self.tfIdTargetId[activeId].append(self.abbrToRnaId[ D["Target"].encode("utf-8") ])
				self.targetIdTFId[self.abbrToRnaId[ D["Target"].encode("utf-8") ]].append(activeId)

			self.tfKd[D["TF"].encode("utf-8")] = D["kd"]
		self.tfNTargets = dict([(key, len(val)) for key,val in self.tfTarget.iteritems()])
		return

	def pPromoterBound(self, dissocConstant, nPromoters, nTfs):
		b = float(dissocConstant + nPromoters + nTfs)
		c = float(nPromoters * nTfs)
		return (b - np.sqrt(b**2 - 4 * c)) / (2 * nPromoters)

	def pTfBound(self, dissocConstant, nPromoters, nTfs):
		b = float(dissocConstant + nPromoters + nTfs)
		c = float(nPromoters * nTfs)
		return (b - np.sqrt(b**2 - 4 * c)) / (2 * nTfs)

	def _buildLookups(self, raw_data, sim_data):
		geneIdToRnaId = dict([(x["geneId"].encode("utf-8"), x["id"].encode("utf-8")) for x in raw_data.rnas])

		self.abbrToRnaId = {}
		for lookupInfo in raw_data.tfIds:
			if len(lookupInfo["geneId"]) == 0:
				continue
			self.abbrToRnaId[lookupInfo["TF"].encode("utf-8")] = geneIdToRnaId[lookupInfo["geneId"].encode("utf-8")]

		self.abbrToActiveId = dict([(x["TF"].encode("utf-8"), x["activeId"].encode("utf-8").split(", ")) for x in raw_data.tfIds if len(x["activeId"]) > 0])

		# Get rid of two-component systems, which we will incorporate later
		keysToDelete = []
		excludes = ["MONOMER0-1",
					"PROTEIN-NRIP",
					"MONOMER0-4180",
					"CPLX0-7754",
					"CPLX0-7884",
					"MONOMER0-4276",
					"CPLX0-7748",
					"MONOMER0-4288",
					"CPLX0-7721",
					"MONOMER0-2741",
					"MONOMER0-4147",
					"CPLX0-7795",
					"MONOMER0-4341",
					]
		for key, val in self.abbrToActiveId.iteritems():
			if val[0].startswith("PHOSPHO") or val[0] in excludes:
				keysToDelete.append(key)
				continue
			locations = sim_data.getter.getLocation(val)
			L = []
			for frameId, loc in zip(val, locations):
				L.append(frameId + "[" + loc[0] + "]")
			self.abbrToActiveId[key] = L

		for key in keysToDelete:
			del self.abbrToActiveId[key]