#!/usr/bin/env python

"""
Replication

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/18/2014
@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/12/2014
"""

from __future__ import division

import numpy as np

import wholecell.processes.process
import wholecell.utils.polymerize_matrix

_NT_ORDER = ['A','T','G','C']
_BASE_PAD_VALUE = ' '

class Replication(wholecell.processes.process.Process):
	""" Replication """

	_name = "Replication"

	# Constructor
	def __init__(self):
		# Constants
		self.sequence = None
		self.genomeLength = None
		self.dnaPolymeraseElongationRate = None

		super(Replication, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(Replication, self).initialize(sim, kb)

		# Load parameters
		dNtpIds = ['DATP[c]', 'DTTP[c]', 'DCTP[c]', 'DGTP[c]']
		dNmpIds = ['DAMP[n]', 'DTMP[n]', 'DCMP[n]', 'DGMP[n]']

		self.sequence = np.array(list(kb.genomeSeq))
		self.genomeLength = kb.genomeLength
		self.dnaPolymeraseElongationRate = kb.dnaPolymeraseElongationRate.to('nucleotide / s').magnitude * self.timeStepSec
		oricCenter = kb.oriCCenter.to('nucleotide').magnitude

		# Views
		self.dntps = self.bulkMoleculesView(dNtpIds)
		self.dnmps = self.bulkMoleculesView(dNmpIds)
		self.ppi = self.bulkMoleculeView('PPI[c]')

		self.h2o = self.bulkMoleculeView('H2O[c]')

		self.dnaPolymerase = self.uniqueMoleculesView('activeDnaPolymerase')
		
		self.dnaPolymerase.moleculeNew('activeDnaPolymerase', chromosomeLocation = oricCenter, directionIsPositive = True)
		self.dnaPolymerase.moleculeNew('activeDnaPolymerase', chromosomeLocation = oricCenter, directionIsPositive = False)

	def calculateRequest(self):
		self.dnaPolymerase.requestAll()

		dnaPolymerase = self.dnaPolymerase.molecules()

		totalNtRequest = np.array([0]*len(_NT_ORDER))
		for p in dnaPolymerase:
			totalNtRequest += self.calculateNucleotideRequest(p)

		# Assumes reaction taking place is:
		# dNTP + H2O --> dNMP + PPi
		self.dntps.requestIs(totalNtRequest)
		self.h2o.requestIs(np.sum(totalNtRequest))

	# Calculate temporal evolution
	def evolveState(self):
		allDnaPolymerase = self.dnaPolymerase.molecules()
		if len(allDnaPolymerase) == 0:
			return

		# Build sequence matrix for polymerize_matrix function
		sequenceMatrix = self.buildSequenceMatrix(allDnaPolymerase)
		ntpCounts = self.dntps.counts()
		bases = np.array(_NT_ORDER)
		energy = 0
		energyCostPerBase = 0
		polymeraseProgress, baseAmounts, baseCosts, energy, energyCost = wholecell.utils.polymerize_matrix.PolymerizeMatrix(
					sequenceMatrix, ntpCounts, bases, _BASE_PAD_VALUE, energy, energyCostPerBase
					)

		# Update DNA polymerase locations based on polymerization polymeraseProgress
		for i,dnaPolymerase in enumerate(allDnaPolymerase):
			self.updatePolymerasePosition(dnaPolymerase, polymeraseProgress[i])
			
		# Update metabolite counts based on polymerization polymeraseProgress
		# Assumes reaction taking place is:
		# dNTP + H2O --> dNMP + PPi
		self.ppi.countInc(np.sum(baseCosts))
		self.dnmps.countsInc(baseCosts)
		self.dntps.countsDec(baseCosts)
		self.h2o.countDec(np.sum(baseCosts))

	def calculateNucleotideRequest(self, dnaPolymerase):
		'''Calculates nucleotide request based on sequence'''
		seq = self.calculateUpcomingSequence(dnaPolymerase)
		return np.array([np.sum(seq == nt) for nt in _NT_ORDER])

	def buildSequenceMatrix(self, allDnaPolymerase):
		'''Builds sequence matrix for polymerize function'''
		sequenceList = []
		for dnaPolymerase in allDnaPolymerase:
			sequenceList.append(self.calculateUpcomingSequence(dnaPolymerase).tolist())
			# TODO: Do this in numpy don't convert to list

		maxLen = max([len(x) for x in sequenceList])

		for s in sequenceList:
			diff = maxLen - len(s)
			if diff > 0:
				s.extend([_BASE_PAD_VALUE]*diff)

		return np.matrix(sequenceList)

	def calculateUpcomingSequence(self, dnaPolymerase):
		'''Wraps actual sequence calculation'''
		return calculateSequence(
			dnaPolymerase.attr('chromosomeLocation'),
			dnaPolymerase.attr('directionIsPositive'),
			self.dnaPolymeraseElongationRate,
			self.sequence,
			self.genomeLength
			)

	def updatePolymerasePosition(self, dnaPolymerase, polymeraseProgress):
		''' Wraps actual update calculation'''
		dnaPolymerase.attrIs(chromosomeLocation = 
					calculatePolymerasePositionUpdate(
						dnaPolymerase.attr('chromosomeLocation'),
						dnaPolymerase.attr('directionIsPositive'),
						polymeraseProgress,
						self.genomeLength
					)
				)

def calculatePolymerasePositionUpdate(currentPosition, directionIsPositive, difference, genomeLength):
	if difference < 0:
		raise Exception, 'Polymerase position difference is negative!\n'

	if directionIsPositive:
		return (currentPosition + difference) % genomeLength
	else:
		return (currentPosition - difference) % genomeLength

def calculateSequence(chromosomeLocation, directionIsPositive, elongationRate, sequence, genomeLength):
	'''
	Calculates sequence in front of DNA polymerase
	based on position, direction, and eloncation rate
	'''

	if directionIsPositive:
		return sequence[
		np.arange(chromosomeLocation, chromosomeLocation + elongationRate) % genomeLength
		]
	else:
		return sequence[
		np.arange(chromosomeLocation, chromosomeLocation - elongationRate, -1) % genomeLength
		]

	# TODO: Check for hitting the end of the chromosome!