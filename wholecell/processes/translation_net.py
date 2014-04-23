#!/usr/bin/env python

"""
TranslationNet

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/23/2014
"""

from __future__ import division

import numpy as np

import wholecell.processes.process

class TranslationNet(wholecell.processes.process.Process):
	""" TranslationNet """

	_name = "TranslationNet"


	# Construct object graph
	def initialize(self, sim, kb):
		super(TranslationNet, self).initialize(sim, kb)

		# Views
		self.aas = self.bulkMoleculesView(aaIDs)


	def calculateRequest(self):
		self.aas.requestAll()


	# Calculate temporal evolution
	def evolveState(self):
		# TODO: implement once we do something other than manipulate AAs
		pass

aaIDs = [
	"ALA-L[c]", "ARG-L[c]", "ASN-L[c]", "ASP-L[c]", "CYS-L[c]", "GLU-L[c]",
	"GLN-L[c]", "GLY[c]", "HIS-L[c]", "ILE-L[c]",  "LEU-L[c]", "LYS-L[c]",
	"MET-L[c]", "PHE-L[c]", "PRO-L[c]", "SELNP[c]", "SER-L[c]", "THR-L[c]",
	"TRP-L[c]", "TYR-L[c]", "VAL-L[c]",
	]
