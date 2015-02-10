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

		self.KcatEndoRNaseFullRNA = kb.KcatEndoRNaseFullRNA.asNumber(1 / units.s) * self.timeStepSec

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
		self.h = self.bulkMoleculesView(['H[c]'])

		self.fragmentMetabolites = self.bulkMoleculesView(endCleavageMetaboliteIds)
		self.fragmentBases = self.bulkMoleculesView([id_ + "[c]" for id_ in kb.fragmentNT_IDs])

		self.endoRnases = self.bulkMoleculesView(endoRnaseIds)
		self.exoRnases = self.bulkMoleculesView(exoRnaseIds)
		self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_DEGRADATION)


	# Calculate temporal evolution

	def calculateRequest(self):
		
		KcatEndoRNaseFragmentMin = 0.005 # cleavages/s
		KcatEndoRNaseFragmentMax = 0.231 # cleavages/s

		RNAspecificity = self.rnaDegRates / self.rnaDegRates.sum()

		nRNAsTotalToDegrade = self.KcatEndoRNaseFullRNA * self.endoRnases.total().sum()

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
		waterForNewRnas = (nRNAsToDegrade * (self.rnaLens - 1)).sum() + nRNAsToDegrade.sum()
		waterForLeftOverFragments = self.fragmentBases.total().sum()
		self.h2o.requestIs(waterForNewRnas + waterForLeftOverFragments)

	def evolveState(self):
		self.writeToListener("RnaDegradationListener", "countRnaDegraded", self.rnas.counts())
		self.writeToListener("RnaDegradationListener", "nucleotidesFromDegradation", (self.rnas.counts() * self.rnaLens).sum())

		KcatEndoRNaseFragmentMin = 0.005 # cleavages/s
		KcatEndoRNaseFragmentMax = 0.231 # cleavages/s

		# Calculate endolytic cleavage events
		# Modeling assumption: Once a RNA is cleaved by an endonuclease it's resulting nucleotides
		# are lumped together as "polymerized fragments". These fragments can carry over from
		# previous timesteps. We are also assuming that during endonucleolytic cleavage the 5'
		# terminal phosphate is removed. This is modeled as all of the fragments being one
		# long linear chain of "fragment bases".
		# Example:
		# PPi-Base-PO4(-)-Base-PO4(-)-Base-PO4(-)-Base-OH
		#			==>
		#			Pi-FragmentBase-PO4(-)-FragmentBase-PO4(-)-FragmentBase-PO4(-)-FragmentBase + PPi
		# Note: Lack of -OH on 3' end of chain
		metabolitesEndoCleavage = np.dot(self.endoDegradationSMatrix, self.rnas.counts())
		self.rnas.countsIs(0)
		self.fragmentMetabolites.countsInc(metabolitesEndoCleavage)

		# Check if can happen exonucleolytic digestion
		if self.fragmentBases.counts().sum() == 0:
			return

		# Calculate exolytic cleavage events
		# Modeling assumption: We model fragments as one long fragment chain of polymerized nucleotides.
		# We are also assuming that there is no sequence specificity or bias towards which nucleotides
		# are hydrolyzed.
		# Example
		# Pi-FragmentBase-PO4(-)-FragmentBase-PO4(-)-FragmentBase-PO4(-)-FragmentBase + 4 H2O
		#			==>
		#			3 NMP + 3 H(+)
		# Note: Lack of -OH on 3' end of chain
		
		kcatExoRNase = 50 # nucleotides/s
		nExoRNases = self.exoRnases.counts()
		exoCapacity = nExoRNases.sum() * kcatExoRNase

		if exoCapacity >= self.fragmentBases.counts().sum():
			self.nmps.countsInc(self.fragmentBases.counts())
			self.h2o.countsDec(self.fragmentBases.counts().sum())
			self.h.countsInc(self.fragmentBases.counts().sum())
			self.fragmentBases.countsIs(0)
		else:
			fragmentSpecificity = self.fragmentBases.counts() / self.fragmentBases.counts().sum()
			possibleBasesToDigest = self.randomState.multinomial(exoCapacity, fragmentSpecificity)
			fragmentBasesDigested = self.fragmentBases.counts() - np.fmax(self.fragmentBases.counts() - possibleBasesToDigest, 0)
			self.nmps.countsInc(fragmentBasesDigested)
			self.h2o.countsDec(fragmentBasesDigested.sum())
			self.h.countsInc(fragmentBasesDigested.sum())
			self.fragmentBases.countsDec(fragmentBasesDigested)
