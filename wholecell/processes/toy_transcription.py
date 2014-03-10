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

		# HACK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

		# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		self.unboundRNApoly = self.uniqueMoleculesView('RNA polymerase',
			boundToChromosome = ('==', False))

		self.boundRNApoly = self.uniqueMoleculesView('RNA polymerase',
			boundToChromosome = ('==', True))


	def calculateRequest(self):
		nToBind = self.unboundRNApoly.total() * self.bindingProb
		nToUnbind = self.boundRNApoly.total() * (1 - self.bindingProb)

		self.unboundRNApoly.requestIs(nToBind)
		self.boundRNApoly.requestIs(nToUnbind)


	# Calculate temporal evolution
	def evolveState(self):
		bound = self.boundRNApoly.molecules()
		unbound = self.unboundRNApoly.molecules()

		for molecule in bound:
			molecule.attrIs('boundToChromosome', False)
			molecule.attrIs('chromosomeLocation', 0)

		for molecule in unbound:
			molecule.attrIs('boundToChromosome', True)
			molecule.attrIs('chromosomeLocation',
				self.randStream.randi(self.chromosomeLength)
				)

		self.unboundRNApoly.moleculesNew(
			'RNA polymerase',
			self.creationRate
			)

		print len(bound), len(unbound)
