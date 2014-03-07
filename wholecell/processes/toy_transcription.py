#!/usr/bin/env python

"""
ToyTranscription

Creates unique RNA polymerases and binds them to an imaginary chromosome.

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/28/2014
"""

from __future__ import division

import numpy as np

import wholecell.processes.process

class ToyTranscription(wholecell.processes.process.Process):
	""" ToyTranscription """

	# Constructor
	def __init__(self):
		self.meta = {
			"id": "ToyTranscription",
			"name": "ToyTranscription",
		}

		self.chromosomeLength = 1000000
		self.initialCounts = 2000
		self.footprint = 20
		self.bindingProb = 0.7
		self.creationRate = 50


		super(ToyTranscription, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(ToyTranscription, self).initialize(sim, kb)

		container = sim.states["UniqueMolecules"]._container

		container.objectsNew('RNA polymerase', self.initialCounts)
		for molecule in container.iterObjects('RNA polymerase'):
			if self.randStream.rand() <= self.bindingProb:
				molecule.attrIs('boundToChromosome', True)
				location = self.randStream.randi(self.chromosomeLength)
				molecule.attrIs('chromosomeLocation', location)

			else:
				molecule.attrIs('boundToChromosome', False)
				molecule.attrIs('chromosomeLocation', -1)

		self.unboundRNApoly_query = sim.states['UniqueMolecules'].queryNew(
			'RNA polymerase', boundToChromosome = ('==', False))

		self.boundRNApoly_query = sim.states['UniqueMolecules'].queryNew(
			'RNA polymerase', boundToChromosome = ('==', True))

	def requestBulkMolecules(self):
		# TODO: Needs to request bulk RNA polymerase molecules
		pass

	def requestUniqueMolecules(self):
		unboundMolecules = self.unboundRNApoly_query.objects()
		boundMolecules = self.boundRNApoly_query.objects()

		nToBind = len(unboundMolecules) * self.bindingProb
		nToUnbind = len(boundMolecules) * (1 - self.bindingProb)

		self.uniqueMoleculesPartition.requestByMolecules(nToBind, unboundMolecules)
		self.uniqueMoleculesPartition.requestByMolecules(nToUnbind, boundMolecules)


	# Calculate temporal evolution
	def evolveState(self):
		# Bind and unbind molecules
		bound = self.uniqueMoleculesPartition.evaluateQuery('RNA polymerase', boundToChromosome = ('==', True))
		unbound = self.uniqueMoleculesPartition.evaluateQuery('RNA polymerase', boundToChromosome = ('==', False))

		for molecule in bound:
			if self.randStream.rand() <= (1. - self.bindingProb):
				molecule.attrIs('boundToChromosome', False)
				molecule.attrIs('chromosomeLocation', -1)

		for molecule in unbound:
			if self.randStream.rand() <= self.bindingProb:
				molecule.attrIs('boundToChromosome', True)
				location = self.randStream.randi(self.chromosomeLength)
				molecule.attrIs('chromosomeLocation', location)

		# Add some new molecules to replenish molecules lost to degradation

		self.uniqueMoleculesPartition.moleculesNew('RNA polymerase', self.creationRate)
