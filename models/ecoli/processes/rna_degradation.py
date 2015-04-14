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
import random

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

		rnaIds = kb.process.transcription.rnaData['id']

		#RNase
		exoRnaseIds = ["EG11620-MONOMER[c]", "G7175-MONOMER[c]", "EG10858-MONOMER[c]",  "EG10863-MONOMER[c]", "EG11259-MONOMER[c]", "EG11547-MONOMER[c]", "EG10746-MONOMER[c]", "EG10743-MONOMER[c]", "G7842-MONOMER[c]"]
		endoRnaseIds = ["EG10856-MONOMER[p]", "EG10857-MONOMER[c]", "G7175-MONOMER[c]", "EG10859-MONOMER[c]", "EG11299-MONOMER[c]", "EG10860-MONOMER[c]", "EG10861-MONOMER[c]", "G7365-MONOMER[c]", "EG10862-MONOMER[c]"]

		self.KcatEndoRNaseFullRNA = kb.constants.KcatEndoRNaseFullRNA.asNumber(1 / units.s) * self.timeStepSec
		# KcatEndoRNasesFullRNA	[0.0015, 0.0015, 0.0015, 0.183, 0.023, 0.023, 0.0015, 0.0015, 0.007]	1/units.s
		self.KcatEndoRNasesFullRNA = kb.constants.KcatEndoRNasesFullRNA.asNumber() * self.timeStepSec

		# Rna
		self.rnaDegRates = kb.process.transcription.rnaData['degRate'].asNumber()
		self.isMRna = kb.process.transcription.rnaData["isMRna"]
		# expectedMRnaDegradationRate = kb.process.transcription.rnaData['degRate'][isMRna].asNumber()
		self.isRRna = kb.process.transcription.rnaData["isRRna"]
		# expectedRRnaDegradationRate = kb.process.transcription.rnaData['degRate'][isRRna].asNumber()
		self.isTRna = kb.process.transcription.rnaData["isTRna"]
		# expectedTRnaDegradationRate = kb.process.transcription.rnaData['degRate'][isTRna].asNumber()

		# import ipdb; ipdb.set_trace()


		self.rnaLens = kb.process.transcription.rnaData['length'].asNumber()

		# Build stoichiometric matrix
		endCleavageMetaboliteIds = [id_ + "[c]" for id_ in kb.moleculeGroups.fragmentNT_IDs]
		endCleavageMetaboliteIds.extend(["H2O[c]", "PPI[c]", "H[c]"])
		nmpIdxs = range(4)
		h2oIdx = endCleavageMetaboliteIds.index("H2O[c]")
		ppiIdx = endCleavageMetaboliteIds.index("PPI[c]")
		hIdx = endCleavageMetaboliteIds.index("H[c]")
		self.endoDegradationSMatrix = np.zeros((len(endCleavageMetaboliteIds), len(rnaIds)), np.int64)
		self.endoDegradationSMatrix[nmpIdxs, :] = units.transpose(kb.process.transcription.rnaData['countsACGU']).asNumber()
		self.endoDegradationSMatrix[h2oIdx, :] = 0
		self.endoDegradationSMatrix[ppiIdx, :] = 1
		self.endoDegradationSMatrix[hIdx, :] = 0
		
		# Views
		self.rnas = self.bulkMoleculesView(rnaIds)
		self.h2o = self.bulkMoleculesView(['H2O[c]'])
		self.nmps = self.bulkMoleculesView(["AMP[c]", "CMP[c]", "GMP[c]", "UMP[c]"])
		self.h = self.bulkMoleculesView(['H[c]'])

		self.fragmentMetabolites = self.bulkMoleculesView(endCleavageMetaboliteIds)
		self.fragmentBases = self.bulkMoleculesView([id_ + "[c]" for id_ in kb.moleculeGroups.fragmentNT_IDs])

		self.endoRnases = self.bulkMoleculesView(endoRnaseIds)
		self.exoRnases = self.bulkMoleculesView(exoRnaseIds)
		self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_DEGRADATION)


	# Calculate temporal evolution

	def calculateRequest(self):

		# Computing RNA specificity according to the measured RNA decays 
		# and accessibility (amount of RNAs)
		TotalDegradationRate = self.rnaDegRates * self.rnas.total() * self.timeStepSec
		# TotalDegradationRate = self.rnaDegRates * self.timeStepSec
		RNAspecificity = TotalDegradationRate / TotalDegradationRate.sum()


		# Calculating fraction of EndoRNases needed 
		print TotalDegradationRate.sum() / self.endoRnases.total().sum()
		print self.endoRnases.total().sum() 
		print TotalDegradationRate.sum()
		FractionActiveEndoRNases = 1
		# import ipdb; ipdb.set_trace()
		# if TotalDegradationRate.sum() < self.KcatEndoRNaseFullRNA * self.endoRnases.total().sum():
		if TotalDegradationRate.sum() < sum(self.KcatEndoRNasesFullRNA * self.endoRnases.total()):
			FractionActiveEndoRNases = TotalDegradationRate.sum() / (
					# self.KcatEndoRNaseFullRNA * self.endoRnases.total().sum()
					sum(self.KcatEndoRNasesFullRNA * self.endoRnases.total())
				)
		print FractionActiveEndoRNases



		# Calculating the total number of RNAs to degrade according to
		# the total number of "active" endoRNases and their cleavage activity
		# nRNAsTotalToDegrade = np.round(self.KcatEndoRNaseFullRNA * 
		nRNAsTotalToDegrade = np.round(sum(self.KcatEndoRNasesFullRNA * 
				# self.endoRnases.total().sum() * 
				self.endoRnases.total()) * 
				FractionActiveEndoRNases
			)
		MrnaSpec = sum(RNAspecificity * self.isMRna)
		TrnaSpec = sum(RNAspecificity * self.isTRna)
		RrnaSpec = sum(RNAspecificity * self.isRRna)
		TargetEndoRNasesFullMRNA = [MrnaSpec, MrnaSpec, 			0, MrnaSpec, MrnaSpec, 				MrnaSpec, MrnaSpec, MrnaSpec, MrnaSpec]
		TargetEndoRNasesFullTRNA = [TrnaSpec, 0, 					1, TrnaSpec, 0, 					TrnaSpec, TrnaSpec, TrnaSpec, TrnaSpec]
		TargetEndoRNasesFullRRNA = [RrnaSpec, TrnaSpec + RrnaSpec, 	0, RrnaSpec, TrnaSpec + RrnaSpec, 	RrnaSpec, RrnaSpec, RrnaSpec, RrnaSpec]
		# import ipdb; ipdb.set_trace()
		nMRNAsTotalToDegrade = np.round(sum(TargetEndoRNasesFullMRNA *
				self.endoRnases.total() * 
				self.KcatEndoRNasesFullRNA) *
				FractionActiveEndoRNases
			)
		nTRNAsTotalToDegrade = np.round(sum(TargetEndoRNasesFullTRNA *
				self.endoRnases.total() * 
				self.KcatEndoRNasesFullRNA) * 
				FractionActiveEndoRNases
			)
		nRRNAsTotalToDegrade = np.round(sum(TargetEndoRNasesFullRRNA *
				self.endoRnases.total() * 
				self.KcatEndoRNasesFullRNA) * 
				FractionActiveEndoRNases
			)
		# print nMRNAsTotalToDegrade + nTRNAsTotalToDegrade + nRRNAsTotalToDegrade
		# print nRNAsTotalToDegrade
		if nRNAsTotalToDegrade != nMRNAsTotalToDegrade + nTRNAsTotalToDegrade + nRRNAsTotalToDegrade:
			nRNAsTotalToDegrade = nMRNAsTotalToDegrade + nTRNAsTotalToDegrade + nRRNAsTotalToDegrade
		# print nRNAsTotalToDegrade


		nRNAsToDegrade = np.zeros(len(RNAspecificity))
		nMRNAsToDegrade = np.zeros(len(RNAspecificity))
		nTRNAsToDegrade = np.zeros(len(RNAspecificity))
		nRRNAsToDegrade = np.zeros(len(RNAspecificity))
		
		# import ipdb; ipdb.set_trace()
		# while nRNAsToDegrade.sum() < nRNAsTotalToDegrade:
			# nRNAsToDegrade += np.fmin(
			# 		self.randomState.multinomial(nRNAsTotalToDegrade - nRNAsToDegrade.sum(), RNAspecificity),
			# 		self.rnas.total()
			# 	)
		nRNAs = self.rnas.total()
		for i in range(0,len(RNAspecificity)):
			if self.rnas.total()[i] != 0:
				nRNAs[i] = 1

		while nMRNAsToDegrade.sum() < nMRNAsTotalToDegrade and (self.rnas.total() * self.isMRna).sum() != 0:
			nMRNAsToDegrade += np.fmin(
					self.randomState.multinomial(nMRNAsTotalToDegrade - nMRNAsToDegrade.sum(), 1. / sum(RNAspecificity * self.isMRna * nRNAs) * RNAspecificity * self.isMRna * nRNAs),
					self.rnas.total() * self.isMRna
				)
		while nTRNAsToDegrade.sum() < nTRNAsTotalToDegrade and (self.rnas.total() * self.isTRna).sum() != 0:
			nTRNAsToDegrade += np.fmin(
					self.randomState.multinomial(nTRNAsTotalToDegrade - nTRNAsToDegrade.sum(), 1. / sum(self.isTRna * nRNAs) * self.isTRna * nRNAs),
					self.rnas.total() * self.isTRna
				)
		while nRRNAsToDegrade.sum() < nRRNAsTotalToDegrade and (self.rnas.total() * self.isRRna).sum() != 0:
			nRRNAsToDegrade += np.fmin(
					self.randomState.multinomial(nRRNAsTotalToDegrade - nRRNAsToDegrade.sum(),  1. / sum(self.isRRna * nRNAs) * self.isRRna * nRNAs),
					self.rnas.total() * self.isRRna
				)

		nRNAsToDegrade = nMRNAsToDegrade + nTRNAsToDegrade + nRRNAsToDegrade
		# print nRNAsToDegrade.sum() / nRNAsTotalToDegrade * 100
		# import ipdb; ipdb.set_trace()
		# print 'real number of RNAs degraded = %f' % nRNAsToDegrade.sum()

		# old method
		# ODE model: dr/dt = kb - kcatEndo * TotalEndoRNases * kd*r/sum_g(kd*r)
		# endoRnasesRelative = self.endoRnases.total().sum() * FractionActiveEndoRNases * RNAspecificity
		# nRNAsToDegrade = np.fmin(
		# 	# self.randomState.poisson(self.KcatEndoRNaseFullRNA * self.timeStepSec * endoRnasesRelative),
		# 	self.randomState.poisson(self.rnaDegRates * self.rnas.total() * self.timeStepSec),
		# 	self.rnas.total()
		# 	)
		# import ipdb; ipdb.set_trace()


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
		self.writeToListener("RnaDegradationListener", "FractionActiveEndoRNases", FractionActiveEndoRNases)

	def evolveState(self):
		print "RNAs degraded = %.3f" % sum(self.rnas.counts())
		for i in range(0,len(self.rnas.counts())):
			if self.rnas.counts()[i] != 0:
				print i

		self.writeToListener("RnaDegradationListener", "countRnaDegraded", self.rnas.counts())
		self.writeToListener("RnaDegradationListener", "nucleotidesFromDegradation", (self.rnas.counts() * self.rnaLens).sum())

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