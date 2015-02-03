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
		# Parameters
		self.rnaLens = None			# RNA lengths
		self.rnaDegRates = None		# RNA degradation rates (1/s)
		self.rnaDegSMat = None		# RNA degradation stoichiometry matrix [metabolite x rna]

		# Views
		self.metabolites = None
		self.rnas = None
		self.rnase = None

		self.Rnases = None
		self.endoRnases = None
		self.exoRnases = None
		self.fragmentBases = None

		self.rnaSequences = None

		super(RnaDegradation, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(RnaDegradation, self).initialize(sim, kb)

		metaboliteIds = ["AMP[c]", "CMP[c]", "GMP[c]", "UMP[c]",
			"H2O[c]", "PPI[c]", "H[c]"]

		nmpIdxs = np.arange(0, 4)
		h2oIdx = metaboliteIds.index('H2O[c]')
		ppiIdx = metaboliteIds.index('PPI[c]')
		hIdx = metaboliteIds.index('H[c]')

		rnaIds = kb.rnaData['id']

		#RNase
		exoRnaseIds = ["EG11620-MONOMER[c]", "G7175-MONOMER[c]", "EG10858-MONOMER[c]",  "EG10863-MONOMER[c]", "EG11259-MONOMER[c]", "EG11547-MONOMER[c]", "EG10746-MONOMER[c]", "EG10743-MONOMER[c]", "G7842-MONOMER[c]"]
		endoRnaseIds = ["EG10856-MONOMER[p]", "EG10857-MONOMER[c]", "G7175-MONOMER[c]", "EG10859-MONOMER[c]", "EG11299-MONOMER[c]", "EG10860-MONOMER[c]", "EG10861-MONOMER[c]", "G7365-MONOMER[c]", "EG10862-MONOMER[c]"]

		# Rna
		self.rnaDegRates = kb.rnaData['degRate'].asNumber()
		self.rnaLens = kb.rnaData['length'].asNumber()

		# Build stoichiometric matrix
		# TODO: account for NTP on 5' end
		self.rnaDegSMat = np.zeros((len(metaboliteIds), len(rnaIds)), np.int64)
		self.rnaDegSMat[nmpIdxs, :] = units.transpose(kb.rnaData['countsACGU']).asNumber()
		# self.rnaDegSMat[h2oIdx, :] = -(self.rnaLens - 1)
		self.rnaDegSMat[h2oIdx, :] = -self.rnaLens # using one additional water to hydrolyze PPI on 5' end
		self.rnaDegSMat[ppiIdx, :] = 1
		self.rnaDegSMat[hIdx, :] = self.rnaLens
		
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

		RNAspecificity = self.rnaDegRates / sum(self.rnaDegRates)

		nEndoRNases = self.endoRnases.total()

		nRNAsTotalToDegrade = KcatEndoRNaseFullRNA * sum(nEndoRNases)
		#nRNAsTotalToDegrade = math.floor(KcatEndoRNaseFullRNA * sum(nEndoRNases)) # Do not use all endoRNases

		# nRNAsToDegrade_old =  np.fmin(
		# 	self.randomState.poisson(self.rnaDegRates * self.rnas.total() * self.timeStepSec),
		# 	self.rnas.total()
		# 	)

		if nRNAsTotalToDegrade == 0:
			return

		nRNAsToDegrade = np.fmin(self.randomState.multinomial(nRNAsTotalToDegrade, 
			RNAspecificity),
			self.rnas.total()
			)

		# import ipdb
		# ipdb.set_trace()
		# nReactions = np.dot(self.rnaLens, nRNAsToDegrade_old)

		# self.h2o.requestIs(nReactions)
		self.rnas.requestIs(nRNAsToDegrade)
		self.rnase.requestAll()
		self.endoRnases.requestAll()
		self.exoRnases.requestAll()
		self.fragmentBases.requestAll()

		print sum(nRNAsToDegrade)
		#print -np.dot(self.rnaDegSMat, nRNAsToDegrade)
		metaboliteUsage = np.fmax(
			-np.dot(self.rnaDegSMat, nRNAsToDegrade),
			0
			)

		#print metaboliteUsage

		# TODO: check H2O, PPI and H usage
		#print metaboliteUsage[4] 

		#print metaboliteUsage
		metaboliteUsage[4] = metaboliteUsage[4] + sum(self.fragmentBases.counts()) # H20
		#print metaboliteUsage[4] 
		#import ipdb
		#ipdb.set_trace()

		print metaboliteUsage
		self.metabolites.requestIs(metaboliteUsage)
		#import ipdb
		#ipdb.set_trace()

	def evolveState(self):

		# Check if RNAse R expressed
		#if self.rnase.count() == 0:
		#	return

		# Degrade RNA
		#self.metabolites.countsInc(np.dot(self.rnaDegSMat, self.rnas.counts()))





		# Add ACGU content of fragments from endonucleolytic cleavages to the previous pull of fragments
		metabolitesEndoCleavage = np.dot(self.rnaDegSMat, self.rnas.counts())
		print sum(self.rnas.counts())
		#print -np.dot(self.rnaDegSMat, self.rnas.counts())

		fragmentACGUCount = metabolitesEndoCleavage[0:4] # TODO according to number of rnas.counts allocated

		fragmentACGUCount = fragmentACGUCount + self.fragmentBases.counts()
		print fragmentACGUCount


		# Check if can happen exonucleolytic digestion
		if sum(fragmentACGUCount) == 0:
			return

		# # Compute exoRNases capacity and fragment specificity
		kcatExoRNase = 50 # nucleotides/s
		nExoRNases = self.exoRnases.counts()
		exoCapacity = sum(nExoRNases) * kcatExoRNase

		fragmentSpecificity = fragmentACGUCount / sum(fragmentACGUCount)


		# # Use exoRNases capacity to degrade fragmentACGUCount
		fragmentExoCapacity = self.randomState.multinomial(exoCapacity, 
			fragmentSpecificity)
		#print fragmentExoCapacity

		fragmentACGUDigested = fragmentACGUCount
		#print fragmentACGUDigested

		fragmentACGUCount = fragmentACGUCount - fragmentExoCapacity

		#print fragmentACGUCount
		flag_digestion = 1;

		for x in range(0, len(fragmentACGUCount)):

			if fragmentACGUCount[x] > 0:
				fragmentACGUDigested[x] = fragmentExoCapacity[x]
				flag_digestion = 0
			
			fragmentACGUCount[x] = 0

		if flag_digestion == 1:
			print 'Full digestion'
		#print fragmentACGUDigested
		#print fragmentACGUCount

		# # Balance of ACGUs generated by endonucleolytic cleavages (new fragments) and exonucleolytic digestion
		self.fragmentBases.countsIs(fragmentACGUCount)


		# # Count hydrolysis co-factors related to the exonucleolytic digestion
		H2Oconsumption = -sum(fragmentACGUDigested)
		PPIconsumption = 0 # check!
		Hproduction = sum(fragmentACGUDigested)
		cofactorHidrolysis = np.array((H2Oconsumption, PPIconsumption, Hproduction))

		
		# # Concatenate ACGU content digested and H2O, PPI and H produced/needed
		balanceMetabolites = numpy.concatenate((fragmentACGUDigested, cofactorHidrolysis))
		
		#print balanceMetabolites
		self.metabolites.countsInc(balanceMetabolites)


		# import ipdb
		# ipdb.set_trace()


		self.writeToListener("RnaDegradationListener", "countRnaDegraded", self.rnas.counts())
		self.writeToListener("RnaDegradationListener", "nucleotidesFromDegradation", (self.rnas.counts() * self.rnaLens).sum())

		self.rnas.countsIs(0) # Is that correct?
