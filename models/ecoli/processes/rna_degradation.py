#!/usr/bin/env python

"""
RnaDegradation

RNA degradation sub-model. Encodes molecular simulation of RNA degradation as two main steps guided by RNases: endonucleolytic cleavage and exonucleolytic digestion 

@author: Javier Carrera
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/26/2015
"""

from __future__ import division

import numpy as np
import math
import numpy

import wholecell.processes.process
from wholecell.utils.constants import REQUEST_PRIORITY_DEGRADATION
from wholecell.utils import units

class RnaDegradation(wholecell.processes.process.Process):
	""" RnaDegradation """

	_name = "RnaDegradation"

	# Constructor
	def __init__(self):
		super(RnaDegradation, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(RnaDegradation, self).initialize(sim, kb)

		rnaIds = kb.rnaData['id']

		#RNase
		exoRnaseIds = ["EG11620-MONOMER[c]", "G7175-MONOMER[c]", "EG10858-MONOMER[c]",  "EG10863-MONOMER[c]", "EG11259-MONOMER[c]", "EG11547-MONOMER[c]", "EG10746-MONOMER[c]", "EG10743-MONOMER[c]", "G7842-MONOMER[c]"]
		endoRnaseIds = ["EG10856-MONOMER[p]", "EG10857-MONOMER[c]", "G7175-MONOMER[c]", "EG10859-MONOMER[c]", "EG11299-MONOMER[c]", "EG10860-MONOMER[c]", "EG10861-MONOMER[c]", "G7365-MONOMER[c]", "EG10862-MONOMER[c]"]

		# Rna
		self.rnaDegRates = kb.rnaData['degRate'].asNumber()
		self.rnaLens = kb.rnaData['length'].asNumber()

		# Build stoichiometric matrix
		endCleavageMetaboliteIds = [id_ + "[c]" for id_ in kb.fragmentNT_IDs]
		endCleavageMetaboliteIds.extend(["H2O[c]", "PPI[c]", "H[c]"])
		nmpIdxs = range(4)
		h2oIdx = endCleavageMetaboliteIds.index("H2O[c]")
		ppiIdx = endCleavageMetaboliteIds.index("PPI[c]")
		hIdx = endCleavageMetaboliteIds.index("H[c]")
		self.endoDegradationSMatrix = np.zeros((len(endCleavageMetaboliteIds), len(rnaIds)), np.int64)
		self.endoDegradationSMatrix[nmpIdxs, :] = units.transpose(kb.rnaData['countsACGU']).asNumber()
		self.endoDegradationSMatrix[h2oIdx, :] = 0
		self.endoDegradationSMatrix[ppiIdx, :] = 1
		self.endoDegradationSMatrix[hIdx, :] = 0
		
		# Views
		self.rnas = self.bulkMoleculesView(rnaIds)
		self.h2o = self.bulkMoleculesView(['H2O[c]'])
		self.nmps = self.bulkMoleculesView(["AMP[c]", "CMP[c]", "GMP[c]", "UMP[c]"])

		self.fragmentMetabolites = self.bulkMoleculesView(endCleavageMetaboliteIds)
		self.fragmentBases = self.bulkMoleculesView([id_ + "[c]" for id_ in kb.fragmentNT_IDs])

		self.endoRnases = self.bulkMoleculesView(endoRnaseIds)
		self.exoRnases = self.bulkMoleculesView(exoRnaseIds)
		self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_DEGRADATION)


	# Calculate temporal evolution

	def calculateRequest(self):
		KcatEndoRNaseFullRNA = 0.002 # cleavages/s
		KcatEndoRNaseFragmentMin = 0.005 # cleavages/s
		KcatEndoRNaseFragmentMax = 0.231 # cleavages/s

		RNAspecificity = self.rnaDegRates / self.rnaDegRates.sum()

		nRNAsTotalToDegrade = KcatEndoRNaseFullRNA * self.endoRnases.total().sum()

		if nRNAsTotalToDegrade == 0:
			return

		nRNAsToDegrade = np.fmin(
			self.randomState.multinomial(nRNAsTotalToDegrade, RNAspecificity),
			self.rnas.total()
			)

		self.rnas.requestIs(nRNAsToDegrade)
		self.endoRnases.requestAll()
		self.exoRnases.requestAll()
		self.fragmentBases.requestAll()

		# Calculating amount of water required for total RNA hydrolysis by endo and
		# exonucleases. Assuming complete hydrolysis for now. One additional water
		# for each RNA to hydrolyze the 5' diphosphate.
		self.h2o.requestIs((nRNAsToDegrade * (self.rnaLens - 1)).sum() + nRNAsToDegrade.sum())

	def evolveState(self):
		self.writeToListener("RnaDegradationListener", "countRnaDegraded", self.rnas.counts())
		self.writeToListener("RnaDegradationListener", "nucleotidesFromDegradation", (self.rnas.counts() * self.rnaLens).sum())

		# Calculate endolytic cleavage events
		# Modeling assumption: Once a RNA is cleaved by an endonuclease it's resulting nucleotides
		# are lumped together as "polymerized fragments". These fragments can carry over from
		# previous timesteps. We are also assuming that during endonucleolytic cleavage the 5'
		# terminal phosphate is removed.

		# Example:
		# Step 1: Hydrolyze the endo bond
		# PPi-Base-PO4(-)-Base-PO4(-)-Base-PO4(-)-Base-OH + H2O
		#			==>
		# 			PPi-Base-PO4(-)-Base-OH + PO4H(-)-Base-PO4(-)-Base-OH
		# Step 2: Remove 5' phosphate of left fragment
		# PPi-Base-PO4(-)-Base-OH + H2O
		#			==>
		#			Pi-Base-PO4(-)-Base-OH + HO4P + H(+)
		# Net reaction:
		# PPi-Base-PO4(-)-Base-PO4(-)-Base-PO4(-)-Base-OH + 2 H2O
		#			==>
		#			Pi-Base-PO4(-)-Base-OH + HO4P + H(+) + PO4H(-)-Base-PO4(-)-Base-OH
		metabolitesEndoCleavage = np.dot(self.endoDegradationSMatrix, self.rnas.counts())
		rnasDegraded = self.rnas.counts().sum()
		self.rnas.countsIs(0)
		self.fragmentMetabolites.countsInc(metabolitesEndoCleavage)

		# Check if can happen exonucleolytic digestion
		if self.fragmentBases.counts().sum() == 0:
			return

		# Calculate exolytic cleavage events
		# Modeling assumption: We are assuming that there are no 5' phosphate groups on
		# fragments. We are also assuming that there is no sequence specificity or bias
		# towards which nucleotides are hydrolyzed.
		# Example
		# PO4H(-)-Base-PO4(-)-Base-PO4(-)-BaseOH + 2 H2O
		#			==>
		#			3 PO4H(-)-Base-OH
		# So in general you need N-1 waters
		
		kcatExoRNase = 50 # nucleotides/s
		nExoRNases = self.exoRnases.counts()
		exoCapacity = nExoRNases.sum() * kcatExoRNase
		if exoCapacity > self.fragmentBases.counts().sum():
			self.nmps.countsInc(self.fragmentBases.counts())
			self.h2o.countsDec(self.fragmentBases.counts().sum() - 1)
			self.fragmentBases.countsIs(0)
			#import ipdb; ipdb.set_trace()
		else:
			print 'HEADS UP!'
			fragmentSpecificity = self.fragmentBases.counts() / self.fragmentBases.counts().sum()
			fragmentBasesDigested = self.randomState.multinomial(exoCapacity, fragmentSpecificity)

