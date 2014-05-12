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
import wholecell.utils.polymerize

_NT_ORDER = ['A','T','G','C']

class Replication(wholecell.processes.process.Process):
	""" Replication """

	_name = "Replication"

	# Constructor
	def __init__(self):
		# Constants
		self.elngRate = None

		super(Replication, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(Replication, self).initialize(sim, kb)

		# Load parameters
		dNtpIds = ['DATP[c]', 'DTTP[c]', 'DCTP[c]', 'DGTP[c]']
		dNmpIds = ['DAMP[n]', 'DTMP[n]', 'DCMP[n]', 'DGMP[n]']

		self.sequence = kb.genomeSeq
		self.genomeLength = len(self.sequence)
		self.dnaPolymeraseElongationRate = 200 # nt/s

		# Views
		self.dntps = self.bulkMoleculesView(dNtpIds)
		self.dnmps = self.bulkMoleculesView(dNmpIds)
		self.ppi = self.bulkMoleculeView('PPI[c]')

		self.h2o = self.bulkMoleculeView('H2O[c]')

		self.dnaPolymerase = self.uniqueMoleculesView('activeDnaPolymerase')
		
		self.dnaPolymerase.moleculeNew('activeDnaPolymerase', chromosomeLocation = 3923882, directionIsPositive = True)
		self.dnaPolymerase.moleculeNew('activeDnaPolymerase', chromosomeLocation = 3923882, directionIsPositive = False)


	def calculateRequest(self):
		self.dnaPolymerase.requestAll()

		dnaPolymerase = self.dnaPolymerase.molecules()

		totalNtRequest = np.array([0]*len(_NT_ORDER))
		for p in dnaPolymerase:
			totalNtRequest += self.calculateNucleotideRequest(p)

		self.dntps.requestIs(totalNtRequest)
		self.h2o.requestIs(np.sum(totalNtRequest))

	# Calculate temporal evolution
	def evolveState(self):
		allDnaPolymerase = self.dnaPolymerase.molecules()
		if len(allDnaPolymerase) == 0:
			return

		sequenceMatrix = self.buildSequenceMatrix(allDnaPolymerase)
		ntpCounts = self.dntps.counts()
		bases = np.array(_NT_ORDER)
		basePadValue = ' '
		energy = 0
		energyCostPerBase = 0
		progress, baseAmounts, baseCosts, energy, energyCost = wholecell.utils.polymerize.Polymerize(
					sequenceMatrix, ntpCounts, bases, basePadValue, energy, energyCostPerBase
					)

		for i,p in enumerate(allDnaPolymerase):
			if p.attr('directionIsPositive') == True:
				p.attrIs(chromosomeLocation = (p.attr('chromosomeLocation') + progress[i]) % self.genomeLength)
			else:
				p.attrIs(chromosomeLocation = (p.attr('chromosomeLocation') - progress[i]) % self.genomeLength)

		self.ppi.countInc(np.sum(baseCosts))

		self.dnmps.countsInc(baseCosts)

		self.dntps.countsDec(baseCosts)

		self.h2o.countDec(np.sum(baseCosts))

	def calculateNucleotideRequest(self, dnaPolymerase):
		seq = self.calculateSequence(dnaPolymerase)
		return np.array([seq.count(nt) for nt in _NT_ORDER])

	def buildSequenceMatrix(self, allDnaPolymerase):
		sequence = []
		for p in allDnaPolymerase:
			sequence.append(list(self.calculateSequence(p)))
		import ipdb; ipdb.set_trace()
		return np.matrix(sequence)

	def calculateSequence(self, dnaPolymerase):
		if dnaPolymerase.attr('directionIsPositive') == True:
			start = dnaPolymerase.attr('chromosomeLocation') - 1
			stop = dnaPolymerase.attr('chromosomeLocation') + self.dnaPolymeraseElongationRate
		else:
			start = dnaPolymerase.attr('chromosomeLocation') - self.dnaPolymeraseElongationRate
			stop = dnaPolymerase.attr('chromosomeLocation')
		return self.sequence[start : stop]

		# TODO: Check for hitting the end of the chromosome!

