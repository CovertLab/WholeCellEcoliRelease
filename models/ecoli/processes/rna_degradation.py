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

		metaboliteIds = ["AMP[c]", "CMP[c]", "GMP[c]", "UMP[c]",
			"H2O[c]", "PI[c]", "H[c]"]

		nmpIdxs = np.arange(0, 4)
		h2oIdx = metaboliteIds.index('H2O[c]')
		piIdx = metaboliteIds.index('PI[c]')
		hIdx = metaboliteIds.index('H[c]')

		rnaIds = kb.rnaData['id']

		#RNase
		exoRnaseIds = ["EG11620-MONOMER[c]", "G7175-MONOMER[c]", "EG10858-MONOMER[c]",  "EG10863-MONOMER[c]", "EG11259-MONOMER[c]", "EG11547-MONOMER[c]", "EG10746-MONOMER[c]", "EG10743-MONOMER[c]", "G7842-MONOMER[c]"]
		endoRnaseIds = ["EG10856-MONOMER[p]", "EG10857-MONOMER[c]", "G7175-MONOMER[c]", "EG10859-MONOMER[c]", "EG11299-MONOMER[c]", "EG10860-MONOMER[c]", "EG10861-MONOMER[c]", "G7365-MONOMER[c]", "EG10862-MONOMER[c]"]

		# Rna
		self.rnaDegRates = kb.rnaData['degRate'].asNumber()
		self.rnaLens = kb.rnaData['length'].asNumber()

		# Build stoichiometric matrix
		self.rnaDegSMat = np.zeros((len(nmpIdxs), len(rnaIds)), np.int64)
		self.rnaDegSMat[nmpIdxs, :] = units.transpose(kb.rnaData['countsACGU']).asNumber()
		#self.rnaDegSMat[h2oIdx, :] = -self.rnaLens # using one additional water to hydrolyze PPI on 5' end
		#self.rnaDegSMat[piIdx, :] = 1
		#self.rnaDegSMat[hIdx, :] = self.rnaLens

		# Build stoichiometric matrix
		self.ExornaDegSMat = self.rnaLens #np.zeros((1, len(rnaIds)), np.int64)
		
		# Views
		self.metabolites = self.bulkMoleculesView(metaboliteIds)		
		self.rnas = self.bulkMoleculesView(rnaIds)
		self.rnase = self.bulkMoleculeView('EG11259-MONOMER[c]')

		self.fragmentBases = self.bulkMoleculesView(['Fragment ADN[c]', 'Fragment CYTD[c]', 'Fragment GSN[c]', 'Fragment URI[c]'])

		self.endoRnases = self.bulkMoleculesView(endoRnaseIds)
		self.exoRnases = self.bulkMoleculesView(exoRnaseIds)

		self.rnaSequences = kb.transcriptionSequences

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
		self.rnase.requestAll()
		self.endoRnases.requestAll()
		self.exoRnases.requestAll()
		self.fragmentBases.requestAll()

		# Calculating amount of water required for total RNA hydrolysis by endo and
		# exonucleases. Assuming complete hydrolysis for now.
		# metaboliteUsage = np.fmax(
		# 	-np.dot(self.rnaDegSMat, nRNAsToDegrade),
		# 	0
		# 	)

		metaboliteUsage = np.zeros(7, np.int64)

		# Complete hydrolysis of RNAs degraded by endunucleases (2 molecucles of H2O required per endo-cleavage)
		metaboliteUsage[4] += 2 * nRNAsTotalToDegrade

		# Complete hydrolysis of RNAs degraded and residual fragments by exonucleases
		# H2O molecules required = total nucleotides - 1
		nucleotideRNAs = np.dot(self.ExornaDegSMat,
			nRNAsToDegrade
			)
		nucleotideFragments = self.fragmentBases.total().sum()

		metaboliteUsage[4] += nucleotideRNAs + nucleotideFragments - 1


		self.metabolites.requestIs(metaboliteUsage)

	def evolveState(self):

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
		#			Base-PO4(-)-Base-OH + 2 HO4P + H(+)
		# Net reaction:
		# PPi-Base-PO4(-)-Base-PO4(-)-Base-PO4(-)-Base-OH + 2 H2O
		#			==>
		#			Pi-Base-PO4(-)-Base-OH + HO4P + H(+) + PO4H(-)-Base-PO4(-)-Base-OH

		metabolitesEndoCleavage = np.dot(self.rnaDegSMat, self.rnas.counts())
		fragmentACGUCount = metabolitesEndoCleavage#[0:4] # TODO according to number of rnas.counts allocated
		fragmentACGUCount = fragmentACGUCount + self.fragmentBases.counts()


		# Check if can happen exonucleolytic digestion
		if fragmentACGUCount.sum() == 0:
			return

		# Calculate exolytic cleavage events
		# Modeling assumption: We are assuming that there are no 5' phosphate groups on
		# fragments. We are also assuming that there is no sequence specificity or bias
		# towards which nucleotides are hydrolyzed.

		# Example:
		# PO4H(-)-Base-PO4(-)-Base-PO4(-)-Base-OH + 2 H2O
		#			==>
		#			3 PO4H(-)-Base-OH
		# So in general you need N-1 waters

		kcatExoRNase = 50 # nucleotides/s
		nExoRNases = self.exoRnases.counts()
		exoCapacity = nExoRNases.sum() * kcatExoRNase

		fragmentSpecificity = fragmentACGUCount / fragmentACGUCount.sum()

		# # Use exoRNases capacity to degrade fragmentACGUCount
		fragmentExoCapacity = self.randomState.multinomial(exoCapacity, 
			fragmentSpecificity)

		fragmentACGUDigested = fragmentACGUCount

		fragmentACGUCount = fragmentACGUCount - fragmentExoCapacity

		for x in range(0, len(fragmentACGUCount)):

			if fragmentACGUCount[x] > 0:
				fragmentACGUDigested[x] = fragmentExoCapacity[x]
			
			fragmentACGUCount[x] = 0

		# # Balance of ACGUs generated by endonucleolytic cleavages (new fragments) and exonucleolytic digestion
		self.fragmentBases.countsIs(fragmentACGUCount)


		# # Count hydrolysis co-factors related to the endocleavages and exolytic digestion
		# Endocleavages
		H2Oconsumption = 2 * self.rnas.counts().sum()
		Hproduction = self.rnas.counts().sum()
		PIproduction = self.rnas.counts().sum()
		# Exolytic digestion
		H2Oconsumption += fragmentACGUDigested.sum() - 1

		cofactorHidrolysis = np.array((-H2Oconsumption, PIproduction, Hproduction))

		
		# # Concatenate ACGU content digested and H2O, PPI and H produced/needed
		balanceMetabolites = numpy.concatenate((fragmentACGUDigested, cofactorHidrolysis))
		
		self.metabolites.countsInc(balanceMetabolites)


		self.writeToListener("RnaDegradationListener", "countRnaDegraded", self.rnas.counts())
		self.writeToListener("RnaDegradationListener", "nucleotidesFromDegradation", (self.rnas.counts() * self.rnaLens).sum())

		self.rnas.countsIs(0)

		#import ipdb; ipdb.set_trace()