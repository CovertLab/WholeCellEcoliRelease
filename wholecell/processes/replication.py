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
import collections

import wholecell.processes.process
from wholecell.utils.polymerize_new import polymerize, PAD_VALUE

# NOTE: the ordering here is take advantage of vectorized operations
NT_SINGLELETTERS = ["A", "C", "G", "T"]
DNTP_IDS = ["DATP[c]", "DCTP[c]", "DGTP[c]", "DTTP[c]"]
DNMP_IDS = ["DAMP[n]", "DCMP[n]", "DGMP[n]", "DTMP[n]"]
N_NT_TYPES = len(NT_SINGLELETTERS)

class Replication(wholecell.processes.process.Process):
	""" Replication """

	_name = "Replication"

	# Constructor
	def __init__(self):
		# Parameters
		self.sequence = None
		self.genomeLength = None
		self.dnaPolymeraseElongationRate = None

		# Views
		self.dntps = None
		self.dnmps = None
		self.ppi = None
		self.h2o = None
		self.dnaPolymerase = None

		super(Replication, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(Replication, self).initialize(sim, kb)

		## Load parameters
		# Create genome sequence out of small integers
		self.sequence = np.empty(len(kb.genomeSeq), np.int8)
		self.ntMapping = collections.OrderedDict([(ntpId, i) for i, ntpId in enumerate(NT_SINGLELETTERS)])
		
		for i,letter in enumerate(kb.genomeSeq):
			self.sequence[i] = self.ntMapping[letter] # Build genome sequence as small integers

		# Load modeling parameters
		self.genomeLength = kb.genomeLength
		self.dnaPolymeraseElongationRate = kb.dnaPolymeraseElongationRate.to('nucleotide / s').magnitude * self.timeStepSec
		self.tercCenter = kb.terCCenter.to('nucleotide').magnitude

		# Load gene data to keep track of copy number
		geneIds = kb.geneData['name']
		self.geneEndCoordinate = kb.geneData['endCoordinate']
		self.bufferedGeneEndCoordinate = np.concatenate(
			[self.geneEndCoordinate - self.genomeLength, self.geneEndCoordinate, self.geneEndCoordinate + self.genomeLength]
			) # Add buffer so indexing with numpy can be taken advantage of

		## Views
		self.dntps = self.bulkMoleculesView(DNTP_IDS)
		self.dnmps = self.bulkMoleculesView(DNMP_IDS)
		self.ppi = self.bulkMoleculeView('PPI[c]')
		self.h2o = self.bulkMoleculeView('H2O[c]')
		
		self.genes = self.bulkChromosomesView(geneIds)

		self.dnaPolymerase = self.uniqueMoleculesView('dnaPolymerase')


	def calculateRequest(self):
		self.dnaPolymerase.requestAll()

		dnaPolymerase = self.dnaPolymerase.molecules()
		
		totalNtRequest = np.array([0]*len(self.ntMapping))
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
		sequences = self.buildSequenceMatrix(allDnaPolymerase)
		dNtpCounts = self.dntps.counts()
		reactionLimit = self.dntps.counts().sum()
		polymeraseProgress, dNtpsUsed, nElongations = polymerize(
			sequences,
			dNtpCounts,
			reactionLimit,
			self.randomState
			)

		# Update DNA polymerase locations based on polymerization polymeraseProgress
		for i,dnaPolymerase in enumerate(allDnaPolymerase):
			self.updatePolymerasePosition(dnaPolymerase, polymeraseProgress[i])
			
		# Update metabolite counts based on polymerization polymeraseProgress
		# Assumes reaction taking place is:
		# dNTP + H2O --> dNMP + PPi
		self.ppi.countInc(np.sum(dNtpsUsed))
		self.dnmps.countsInc(dNtpsUsed)
		self.dntps.countsDec(dNtpsUsed)
		self.h2o.countDec(np.sum(dNtpsUsed))


	def calculateNucleotideRequest(self, dnaPolymerase):
		'''Calculates nucleotide request based on sequence'''
		seq = self.calculateUpcomingSequence(dnaPolymerase)
		return np.bincount(seq, minlength = N_NT_TYPES)


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
				s.extend([PAD_VALUE]*diff)

		return np.array(sequenceList, dtype=np.int8)


	def calculateUpcomingSequence(self, dnaPolymerase):
		'''Wraps actual sequence calculation'''

		leadingSequence = calculateSequence(
				dnaPolymerase.attr('chromosomeLocation'),
				dnaPolymerase.attr('directionIsPositive'),
				self.dnaPolymeraseElongationRate,
				self.sequence,
				self.genomeLength
				)

		if dnaPolymerase.attr('isLeading'):
			return leadingSequence

		elif not dnaPolymerase.attr('isLeading'):
			return self.reverseComplement(leadingSequence)


	def updatePolymerasePosition(self, dnaPolymerase, polymeraseProgress):
		''' Wraps actual update calculation'''

		self.updateGeneCopynumber(
			dnaPolymerase.attr('chromosomeLocation'),
			dnaPolymerase.attr('directionIsPositive'),
			polymeraseProgress
			)

		dnaPolymerase.attrIs(chromosomeLocation = 
					calculatePolymerasePositionUpdate(
						dnaPolymerase.attr('chromosomeLocation'),
						dnaPolymerase.attr('directionIsPositive'),
						polymeraseProgress,
						self.genomeLength
					)
				)


	def updateGeneCopynumber(self, currentPosition, directionIsPositive, difference):
		'''
		Returns indicies of genes replicated by polymerase based on position and progress of polymerization
		'''

		if directionIsPositive:
			finalLocation = currentPosition + difference
			bufferedReplicatedGenes = (
				self.bufferedGeneEndCoordinate > currentPosition) & (self.bufferedGeneEndCoordinate <= finalLocation
				)

		else:
			finalLocation = currentPosition - difference
			bufferedReplicatedGenes = (
				self.bufferedGeneEndCoordinate >= finalLocation) & (self.bufferedGeneEndCoordinate < currentPosition
				)

		actualReplicatedGenes = bufferedReplicatedGenes.reshape(3,-1).any(0)
		self.genes.countsInc(actualReplicatedGenes)


	def reverseComplement(self, sequenceVector):
		return (N_NT_TYPES - 1) - sequenceVector


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