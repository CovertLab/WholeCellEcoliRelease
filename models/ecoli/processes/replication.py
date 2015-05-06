#!/usr/bin/env python

"""
Replication

TODO:
- build actual chromosome molecules

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
from wholecell.utils.polymerize import polymerize, PAD_VALUE
from wholecell.utils import units

# NOTE: the ordering here is take advantage of vectorized operations
NT_SINGLELETTERS = ["A", "C", "G", "T"]
N_NT_TYPES = len(NT_SINGLELETTERS)

class Replication(wholecell.processes.process.Process):
	""" Replication """

	_name = "Replication"

	# Constructor
	def __init__(self):
		# Parameters
		self.genomeSequence = None
		self.genomeLength = None
		self.dnaPolymeraseElongationRate = None
		self.ntMapping = None
		self.tercCenter = None
		self.geneEndCoordinate = None
		self.bufferedGeneEndCoordinate = None

		# Views
		self.dntps = None
		self.polymerized = None
		self.ppi = None
		self.dnaPolymerase = None

		super(Replication, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(Replication, self).initialize(sim, kb)

		## Load parameters
		# Create genome sequence out of small integers
		self.genomeSequence = np.empty(len(kb.process.replication.genome_sequence), np.int8)
		self.ntMapping = collections.OrderedDict([(ntpId, i) for i, ntpId in enumerate(NT_SINGLELETTERS)])
		
		for i,letter in enumerate(kb.process.replication.genome_sequence):
			self.genomeSequence[i] = self.ntMapping[letter] # Build genome sequence as small integers

		# Load modeling parameters
		self.genomeLength = len(kb.process.replication.genome_sequence)
		self.dnaPolymeraseElongationRate = kb.constants.dnaPolymeraseElongationRate.asNumber(units.nt / units.s) * self.timeStepSec
		self.tercCenter = kb.constants.terCCenter.asNumber(units.nt)
		self.n_completed_dnaPolymerases = 0

		# Load gene data to keep track of copy number
		geneIds = kb.process.replication.geneData['name']
		self.geneEndCoordinate = kb.process.replication.geneData['endCoordinate']
		self.bufferedGeneEndCoordinate = np.concatenate(
			[self.geneEndCoordinate - self.genomeLength, self.geneEndCoordinate, self.geneEndCoordinate + self.genomeLength]
			) # Add buffer so indexing with numpy can be taken advantage of

		## Views
		self.dntps = self.bulkMoleculesView(kb.moleculeGroups.dNtpIds)
		self.polymerized = self.bulkMoleculesView([id_ + "[c]" for id_ in kb.moleculeGroups.polymerizedDNT_IDs])
		self.ppi = self.bulkMoleculeView('PPI[c]')
		
		self.genes = self.bulkChromosomesView(geneIds)

		self.dnaPolymerase = self.uniqueMoleculesView('dnaPolymerase')


	def calculateRequest(self):
		# Request all dna polymerases
		self.dnaPolymerase.requestAll()

		# Look at all current polymerase properties and calcualate
		# dNTP demands based on upcoming sequence
		allDnaPolymerase = self.dnaPolymerase.allMolecules()
		if len(allDnaPolymerase) == 0:
			return
		nPolymerase = len(allDnaPolymerase)
		allChromosomeLocation, allDirectionIsPositive, allIsLeading = allDnaPolymerase.attrs(
			"chromosomeLocation",
			"directionIsPositive",
			"isLeading"
			)
		
		totalNtRequest = np.array([0]*len(self.ntMapping))
		for i in range(nPolymerase):
			totalNtRequest += np.bincount(
				calculateUpcomingSequence(
					allChromosomeLocation[i],
					allDirectionIsPositive[i],
					allIsLeading[i],
					self.dnaPolymeraseElongationRate,
					self.genomeLength,
					self.genomeSequence,
					self.tercCenter
					), minlength = N_NT_TYPES
				)

		# Assumes reaction taking place is:
		# dNTP + H2O --> dNMP + PPi + H
		self.dntps.requestIs(totalNtRequest)

	# Calculate temporal evolution
	def evolveState(self):
		# Get polymerase properties
		for dnaPolymerase in self.dnaPolymerase.molecules():
			if dnaPolymerase.attr('chromosomeLocation') == self.tercCenter:
				self.dnaPolymerase.moleculeDel(dnaPolymerase)
				self.n_completed_dnaPolymerases += 1

		allDnaPolymerase = self.dnaPolymerase.molecules()
		if len(allDnaPolymerase) == 0:
			# TODO: This assumes we only have one round of replication going on!
			if self.n_completed_dnaPolymerases == 4:
				self._sim.dnaReplicationComplete()
			return

		nPolymerase = len(allDnaPolymerase)
		allChromosomeLocation, allDirectionIsPositive, allIsLeading = allDnaPolymerase.attrs(
			"chromosomeLocation",
			"directionIsPositive",
			"isLeading"
			)

		# Build sequence matrix for polymerize_matrix function
		sequences = buildSequenceMatrix(
			nPolymerase,
			allChromosomeLocation,
			allDirectionIsPositive,
			allIsLeading,
			self.dnaPolymeraseElongationRate,
			self.genomeLength,
			self.genomeSequence,
			self.tercCenter
			)

		# Calculate polymerase progress
		dNtpCounts = self.dntps.counts()
		reactionLimit = self.dntps.counts().sum()
		polymeraseProgress, dNtpsUsed, nElongations = polymerize(
			sequences,
			dNtpCounts,
			reactionLimit,
			self.randomState
			)

		# Update DNA polymerase locations based on polymerization polymeraseProgress
		# Update gene copy number based on leading strand polymerase position
		for i in range(nPolymerase):
			if allIsLeading[i]:
				replicatedGenes = calculateReplicatedGenes(
					allChromosomeLocation[i],
					allDirectionIsPositive[i],
					polymeraseProgress[i],
					self.bufferedGeneEndCoordinate
					)
				if replicatedGenes.any():
					self.genes.countsInc(replicatedGenes)

			allDnaPolymerase[i].attrIs(chromosomeLocation = 
						calculatePolymerasePositionUpdate(
							allChromosomeLocation[i],
							allDirectionIsPositive[i],
							polymeraseProgress[i],
							self.genomeLength
						)
					)
				
		# Update metabolite counts based on polymerization polymeraseProgress
		# Assumes reaction taking place is:
		# dNTP + H2O --> dNMP + PPi + H
		self.ppi.countInc(np.sum(dNtpsUsed))
		self.polymerized.countsInc(dNtpsUsed)
		self.dntps.countsDec(dNtpsUsed)

def buildSequenceMatrix(nPolymerase, allChromosomeLocation, allDirectionIsPositive, allIsLeading, dnaPolymeraseElongationRate, genomeLength, genomeSequence, tercCenter):
	'''Builds sequence matrix for polymerize function'''
	sequenceMatrix = np.empty((nPolymerase, dnaPolymeraseElongationRate), np.int8)
	sequenceMatrix.fill(PAD_VALUE)

	for dnaPolyIndex in range(nPolymerase):
		upcomingSequence = calculateUpcomingSequence(
			allChromosomeLocation[dnaPolyIndex],
			allDirectionIsPositive[dnaPolyIndex],
			allIsLeading[dnaPolyIndex],
			dnaPolymeraseElongationRate,
			genomeLength,
			genomeSequence,
			tercCenter
			)
		sequenceMatrix[dnaPolyIndex, :upcomingSequence.shape[0]] = upcomingSequence

	return sequenceMatrix

def calculateUpcomingSequence(chromosomeLocation, directionIsPositive, isLeading, elongationRate, genomeLength, genomeSequence, tercCenter):
	'''Wraps actual sequence calculation'''

	elongationLength = elongationRate

	# Calculate if terC is in next polymerization step and stop polymerase there
	if directionIsPositive:
		upcomingPositions = np.arange(chromosomeLocation, chromosomeLocation + elongationRate) % genomeLength
	else:
		upcomingPositions = np.arange(chromosomeLocation, chromosomeLocation - elongationRate, -1) % genomeLength
		
	if tercCenter in upcomingPositions:
		elongationLength = np.where(upcomingPositions == tercCenter)[0][0]

	if directionIsPositive:
		leadingSequence = genomeSequence[
			np.arange(chromosomeLocation, chromosomeLocation + elongationLength) % genomeLength
			]
	else:
		leadingSequence = genomeSequence[
			np.arange(chromosomeLocation, chromosomeLocation - elongationLength, -1) % genomeLength
			]

	# Return coding or non-coding strand
	if isLeading:
		return leadingSequence
	else:
		return reverseComplement(leadingSequence)
	
def calculateReplicatedGenes(currentPosition, directionIsPositive, difference, bufferedGeneEndCoordinate):
	'''
	Returns indicies of genes replicated by polymerase based on position and progress of polymerization
	'''

	if directionIsPositive:
		finalLocation = currentPosition + difference
		bufferedReplicatedGenes = (
			bufferedGeneEndCoordinate > currentPosition) & (bufferedGeneEndCoordinate <= finalLocation
			)

	else:
		finalLocation = currentPosition - difference
		bufferedReplicatedGenes = (
			bufferedGeneEndCoordinate >= finalLocation) & (bufferedGeneEndCoordinate < currentPosition
			)

	actualReplicatedGenes = bufferedReplicatedGenes.reshape(3,-1).any(0)

	return actualReplicatedGenes


def reverseComplement(sequenceVector):
	return (N_NT_TYPES - 1) - sequenceVector


def calculatePolymerasePositionUpdate(currentPosition, directionIsPositive, difference, genomeLength):
	if difference < 0:
		raise Exception, 'Polymerase position difference is negative!\n'

	if directionIsPositive:
		return (currentPosition + difference) % genomeLength
		
	else:
		return (currentPosition - difference) % genomeLength