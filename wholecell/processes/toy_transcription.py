#!/usr/bin/env python

"""
ToyTranscription

Creates unique RNA polymerases and binds them to an imaginary chromosome.

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/28/2014
"""

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

		self.container = None


		super(ToyTranscription, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(ToyTranscription, self).initialize(sim, kb)

		uniqueMolecules = sim.states["UniqueMolecules"]
		self.container = uniqueMolecules._container

		self.container.moleculesNew('RNA polymerase', self.initialCounts)
		for molecule in self.container.iterMolecules('RNA polymerase'):
			if self.randStream.rand <= self.bindingProb:
				molecule.attrIs('boundToChromosome', True)
				location = self.randStream.randi(self.chromosomeLength)
				molecule.attrIs('chromosomeLocation', location)
			else:
				molecule.attrIs('boundToChromosome', False)
				molecule.attrIs('chromosomeLocation', -1)

	def requestMoleculeCounts(self):
		# TODO: Needs to request bulk RNA polymerase molecules
		pass


	# Calculate temporal evolution
	def evolveState(self):
		bound = self.container.evaluateQuery('RNA polymerase', boundToChromosome = ('==', True))
		unbound = self.container.evaluateQuery('RNA polymerase', boundToChromosome = ('==', False))

		for molecule in bound:
			if self.randStream.rand <= (1. - self.bindingProb):
				molecule.attrIs('boundToChromosome', False)
				molecule.attrIs('chromosomeLocation', -1)

		for molecule in unbound:
			if self.randStream.rand <= self.bindingProb:
				molecule.attrIs('boundToChromosome', True)
				location = self.randStream.randi(self.chromosomeLength)
				molecule.attrIs('chromosomeLocation', location)
