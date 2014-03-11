#!/usr/bin/env python

"""
ToyReplication

A toy process for developing interactions between processes and the Chromosome
State.

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/11/2014
"""

from __future__ import division

import numpy as np

import wholecell.processes.process

class ToyReplication(wholecell.processes.process.Process):
	""" ToyReplication """

	# Constructor
	def __init__(self):
		self.meta = {
			"id": "ToyReplication",
			"name": "ToyReplication",
			}

		self.dnaPolymeraseFootprint = 50
		self.elongationRate = 100

		super(ToyReplication, self).__init__()


	# Construct object graph
	def initialize(self, sim, kb):
		super(ToyReplication, self).initialize(sim, kb)

		# self.chromosomeView = self.chromosomeView()

		self.replicationForks = self.replicationForkView(
			extentAlongParent = self.dnaPolymeraseFootprint + self.elongationRate,
			extentAlongChildren = self.elongationRate
			)


	def calculateRequest(self):
		self.replicationForks.requestAll()


	# Calculate temporal evolution
	def evolveState(self):
		# can't loop over regions, may overlap
		# need a way to handle replication termination
		# probably won't be able to avoid working with numbers for locations on the chromosome
		for region in self.replicationForks.regions():
			(dnaPolymerase,) = region.molecules(moleculeName = 'DNA polymerase') # singleton unpacking
			(forkLocation,) = region.forkLocations() # singleton unpacking

			moveDistance = region.extentAlongParent()

			newDnaPolyLocation = region.moleculeLocation(dnaPolymerase) + moveDistance - self.dnaPolymeraseFootprint

			region.moleculeMove(dnaPolymerase, newDnaPolyLocation)

			region.forkMove(forkLocation, forkLocation + moveDistance)


		dnaPolymerases = self.replicationForks.molecules(moleculeName = 'DNA polymerase')

		forkLocations = self.replicationForks.forkLocations()

		# how to find forks? associate with polymerases?
		# - find fork, then find nearby molecules
		# how to move forks?

		# communicating distance that can be moved
		# origin for molecule footprints

