#!/usr/bin/env python

"""
RnaDegradation
RNA degradation sub-model. 

Mathematical formulation:

dr/dt = Kb - Kd * r
or,

dr/dt = Kb - kcatEndoRNase * EndoRNase * r / (Km + r)
or,

dr/dt = Kb - kcatEndoRNase * EndoRNase * r/Km / (1 + Sum(r/Km))

	where	r = RNA counts
			Kb = RNA production given a RNAP synthesis rate 
			tau = doubling time
			kcatEndoRNase = enzymatic activity for EndoRNases
			kd = RNA degradation rates 
			Km = Michaelis-Menten constants fitted to recapitulate first-order RNA decay:
				kd * r = kcatEndoRNase * EndoRNase * r / (Km + r), non-cooperative EndoRNases
				kd * r = kcatEndoRNase * EndoRNase * r/Km / (1 + sum(r/Km)), cooperation

This sub-model encodes molecular simulation of RNA degradation as two main steps guided by RNases, "endonucleolytic cleavage" and "exonucleolytic digestion":
1. Compute total counts of RNA to be degraded (D) and total capacity for endo-cleavage (C) at each time point
2. Evaluate C and D. If C > D, then define a fraction of active endoRNases 
3. Dissect RNA degraded into different species (mRNA, tRNA, and rRNA) by accounting endoRNases specificity
4. Update RNA fragments (assumption: fragments are represented as a pull of nucleotides) because of endonucleolytic cleavage
5. Compute total capacity of exoRNases and determine fraction of nucleotides that can be diggested
6. Update pull of metabolites (H and H2O) because of exonucleolytic digestion

@author: Javier Carrera
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/26/2015 - Updated 8/10/2015
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
	def initialize(self, sim, sim_data):
		super(RnaDegradation, self).initialize(sim, sim_data)

		rnaIds = sim_data.process.transcription.rnaData['id']

		# Load constants
		self.nAvogadro = sim_data.constants.nAvogadro
		self.cellDensity = sim_data.constants.cellDensity

		#RNase
		endoRnaseIds = sim_data.process.rna_decay.endoRnaseIds
		exoRnaseIds = sim_data.moleculeGroups.exoRnaseIds
		self.KcatExoRNase = sim_data.constants.KcatExoRNase
		self.KcatEndoRNases = sim_data.process.rna_decay.kcats

		self.TargetEndoRNasesFullMRNA_indexes = sim_data.process.rna_decay.TargetEndoRNasesFullMRNA
		self.TargetEndoRNasesFullTRNA_indexes = sim_data.process.rna_decay.TargetEndoRNasesFullTRNA
		self.TargetEndoRNasesFullRRNA_indexes = sim_data.process.rna_decay.TargetEndoRNasesFullRRNA

		self.mrna_index = sim_data.process.rna_decay.mrna_index
		self.trna_index = sim_data.process.rna_decay.trna_index
		self.rrna_index = sim_data.process.rna_decay.rrna_index
		self.rtrna_index = sim_data.process.rna_decay.rtrna_index

		# Rna
		self.rnaDegRates = sim_data.process.transcription.rnaData['degRate']
		self.isMRna = sim_data.process.transcription.rnaData["isMRna"]
		self.isRRna = sim_data.process.transcription.rnaData["isRRna"]
		self.isTRna = sim_data.process.transcription.rnaData["isTRna"]

		self.rnaLens = sim_data.process.transcription.rnaData['length'].asNumber()

		# Build stoichiometric matrix
		endCleavageMetaboliteIds = [id_ + "[c]" for id_ in sim_data.moleculeGroups.fragmentNT_IDs]
		endCleavageMetaboliteIds.extend(["WATER[c]", "PPI[c]", "PROTON[c]"])
		nmpIdxs = range(4)
		h2oIdx = endCleavageMetaboliteIds.index("WATER[c]")
		ppiIdx = endCleavageMetaboliteIds.index("PPI[c]")
		hIdx = endCleavageMetaboliteIds.index("PROTON[c]")
		self.endoDegradationSMatrix = np.zeros((len(endCleavageMetaboliteIds), len(rnaIds)), np.int64)
		self.endoDegradationSMatrix[nmpIdxs, :] = units.transpose(sim_data.process.transcription.rnaData['countsACGU']).asNumber()
		self.endoDegradationSMatrix[h2oIdx, :] = 0
		self.endoDegradationSMatrix[ppiIdx, :] = 1
		self.endoDegradationSMatrix[hIdx, :] = 0
		
		# Views
		self.rnas = self.bulkMoleculesView(rnaIds)
		self.h2o = self.bulkMoleculesView(["WATER[c]"])
		self.nmps = self.bulkMoleculesView(["AMP[c]", "CMP[c]", "GMP[c]", "UMP[c]"])
		self.h = self.bulkMoleculesView(["PROTON[c]"])

		self.fragmentMetabolites = self.bulkMoleculesView(endCleavageMetaboliteIds)
		self.fragmentBases = self.bulkMoleculesView([id_ + "[c]" for id_ in sim_data.moleculeGroups.fragmentNT_IDs])

		self.endoRnases = self.bulkMoleculesView(endoRnaseIds)
		self.exoRnases = self.bulkMoleculesView(exoRnaseIds)
		self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_DEGRADATION)

		self.Km = sim_data.process.transcription.rnaData["KmEndoRNase"]
		self.EndoRNaseCoop = sim_data.constants.EndoRNaseCooperation
		self.EndoRNaseFunc = sim_data.constants.EndoRNaseFunction


	# Calculate temporal evolution

	def calculateRequest(self):

		# load constants 
		cellMass = (self.readFromListener("Mass", "cellMass") * units.fg)
		cellVolume = cellMass / self.cellDensity
		countsToMolar = 1 / (self.nAvogadro * cellVolume)

		# fraction saturated based on Michaelis-Menten kinetics
		if not self.EndoRNaseCoop:
			fracEndoRnaseSaturated = countsToMolar * self.rnas.total() / (self.Km + (countsToMolar * self.rnas.total()))

		# fraction saturated based on generalized Michaelis-Menten kinetics (EndoRNase cooperation)
		if self.EndoRNaseCoop:
			fracEndoRnaseSaturated = (countsToMolar * self.rnas.total()) / self.Km / (1 + units.sum((countsToMolar * self.rnas.total()) / self.Km))
		
		Kd = self.rnaDegRates
		Kcat = self.KcatEndoRNases
		EndoR = sum(self.endoRnases.total())
		RNA = self.rnas.total()
		FractDiffRNAdecay = units.sum(units.abs(Kd * RNA - units.sum(self.KcatEndoRNases * self.endoRnases.total()) * fracEndoRnaseSaturated))
		FractEndoRRnaCounts = EndoR.astype(float) / sum(RNA.astype(float))

		# Calculate total counts of RNAs to degrade according to
		# the total counts of "active" endoRNases and their cleavage activity
		nRNAsTotalToDegrade = np.round((units.sum(self.KcatEndoRNases * self.endoRnases.total()) * fracEndoRnaseSaturated * (units.s * self.timeStepSec())).asNumber().sum())
		
		# Dissect RNA specificity into mRNA, tRNA, and rRNA as well as specific RNases
		MrnaSpec = units.sum(fracEndoRnaseSaturated * self.isMRna)
		TrnaSpec = units.sum(fracEndoRnaseSaturated * self.isTRna)
		RrnaSpec = units.sum(fracEndoRnaseSaturated * self.isRRna)

		TargetEndoRNasesFullMRNA = np.zeros_like(self.TargetEndoRNasesFullMRNA_indexes)
		TargetEndoRNasesFullMRNA[self.TargetEndoRNasesFullMRNA_indexes == self.mrna_index] = MrnaSpec
		TargetEndoRNasesFullMRNA[self.TargetEndoRNasesFullMRNA_indexes == 0] = 0
		TargetEndoRNasesFullMRNA[self.TargetEndoRNasesFullMRNA_indexes == 1] = 1

		TargetEndoRNasesFullTRNA = np.zeros_like(self.TargetEndoRNasesFullTRNA_indexes)
		TargetEndoRNasesFullTRNA[self.TargetEndoRNasesFullTRNA_indexes == self.trna_index] = TrnaSpec
		TargetEndoRNasesFullTRNA[self.TargetEndoRNasesFullTRNA_indexes == 0] = 0
		TargetEndoRNasesFullTRNA[self.TargetEndoRNasesFullTRNA_indexes == 1] = 1

		TargetEndoRNasesFullRRNA = np.zeros_like(self.TargetEndoRNasesFullRRNA_indexes)
		TargetEndoRNasesFullRRNA[self.TargetEndoRNasesFullRRNA_indexes == self.rrna_index] = RrnaSpec
		TargetEndoRNasesFullRRNA[self.TargetEndoRNasesFullRRNA_indexes == self.rtrna_index] = RrnaSpec + TrnaSpec
		TargetEndoRNasesFullRRNA[self.TargetEndoRNasesFullRRNA_indexes == 0] = 0
		TargetEndoRNasesFullRRNA[self.TargetEndoRNasesFullRRNA_indexes == 1] = 1


		# Dissect total counts of RNA degraded into mRNA, tRNA, and rRNA 
		TargetEndoRNasesFullMRNA = MrnaSpec
		TargetEndoRNasesFullTRNA = TrnaSpec
		TargetEndoRNasesFullRRNA = RrnaSpec

		nMRNAsTotalToDegrade = np.round(sum(TargetEndoRNasesFullMRNA *
				self.endoRnases.total() * 
				self.KcatEndoRNases * (units.s * self.timeStepSec())).asNumber()
			)
		nTRNAsTotalToDegrade = np.round(sum(TargetEndoRNasesFullTRNA *
				self.endoRnases.total() * 
				self.KcatEndoRNases * (units.s * self.timeStepSec())).asNumber()
			)
		nRRNAsTotalToDegrade = np.round(sum(TargetEndoRNasesFullRRNA *
				self.endoRnases.total() * 
				self.KcatEndoRNases * (units.s * self.timeStepSec())).asNumber()
			)

		if nRNAsTotalToDegrade != nMRNAsTotalToDegrade + nTRNAsTotalToDegrade + nRRNAsTotalToDegrade:
			nRNAsTotalToDegrade = nMRNAsTotalToDegrade + nTRNAsTotalToDegrade + nRRNAsTotalToDegrade

		# define RNA specificity across genes
		RNAspecificity = (fracEndoRnaseSaturated / units.sum(fracEndoRnaseSaturated)).asNumber()

		nRNAsToDegrade = np.zeros(len(RNAspecificity))
		nMRNAsToDegrade = np.zeros(len(RNAspecificity))
		nTRNAsToDegrade = np.zeros(len(RNAspecificity))
		nRRNAsToDegrade = np.zeros(len(RNAspecificity))
		
		# boolean variable (nRNAs) to track availability of RNAs for a given gene
		nRNAs = self.rnas.total().astype(np.bool)


		# determine mRNAs to be degraded according to RNA specificities and total counts of mRNAs degraded
		while nMRNAsToDegrade.sum() < nMRNAsTotalToDegrade and (self.rnas.total() * self.isMRna).sum() != 0:
			nMRNAsToDegrade += np.fmin(
					self.randomState.multinomial(nMRNAsTotalToDegrade - nMRNAsToDegrade.sum(), 1. / sum(RNAspecificity * self.isMRna * nRNAs) * RNAspecificity * self.isMRna * nRNAs),
					self.rnas.total() * self.isMRna
				)

		# determine tRNAs and rRNAs to be degraded (with equal specificity) depending on total counts degraded, respectively 
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

		# First order decay with non-functional EndoRNase activity 
		# Determine mRNAs to be degraded according to Poisson distribution (Kdeg * RNA)
		if not self.EndoRNaseFunc:
			nRNAsToDegrade = np.fmin(
				self.randomState.poisson( (self.rnaDegRates * self.rnas.total()).asNumber() ),
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

		self.writeToListener("RnaDegradationListener", "FractionActiveEndoRNases", sum(fracEndoRnaseSaturated))
		self.writeToListener("RnaDegradationListener", "DiffRelativeFirstOrderDecay", FractDiffRNAdecay.asNumber())
		self.writeToListener("RnaDegradationListener", "FractEndoRRnaCounts", FractEndoRRnaCounts)

	def evolveState(self):

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
		nExoRNases = self.exoRnases.counts()
		exoCapacity = nExoRNases.sum() * self.KcatExoRNase * (units.s * self.timeStepSec())

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