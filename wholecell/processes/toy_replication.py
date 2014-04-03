#!/usr/bin/env python

"""
ToyReplication

A toy process for developing interactions between processes and the Chromosome
State.

@author: John Mason
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

		self.dnaPolyForwardFootprint = 50
		self.dnaPolyReverseFootprint = 50
		self.elongationRate = 1000

		super(ToyReplication, self).__init__()


	# Construct object graph
	def initialize(self, sim, kb):
		super(ToyReplication, self).initialize(sim, kb)

		self.replicationForksRegions = self.chromosomeForksView( # special view for forks
			extentForward = self.elongationRate + self.dnaPolyForwardFootprint,
			extentReverse = self.dnaPolyReverseFootprint,
			includeMoleculesOnEnds = True
			)


	def calculateRequest(self):
		self.replicationForksRegions.requestAll()


	# Calculate temporal evolution
	def evolveState(self):
		# TODO: check whether any forks would collide and merge
		# TODO: update fork positions in sync (polymerize function)

		# Iterate over each region set
		for region in self.replicationForksRegions.parentRegions():
			# Get the fork in the region
			forks = list(self.replicationForksRegions.forksInRegion(region))

			if len(forks) != 1:
				print 'Invalid number of forks'
				continue

			fork = forks.pop()

			# Get the molecule on the fork
			dnaPoly = self.replicationForksRegions.moleculeBoundOnFork(fork)

			if dnaPoly is None or dnaPoly.name() != 'DNA polymerase':
				print 'No DNA polymerase found on fork'
				continue

			# Find the furthest extent that the fork can be moved
			extent = self.replicationForksRegions.maximumExtentPastFork(
				fork = fork,
				maxExtentForward = self.elongationRate + self.dnaPolyForwardFootprint,
				)

			# Find any non-DNA poly molecules in the extent
			nonPolymeraseMolecules = set(self.replicationForksRegions.moleculesBoundPastFork(
				fork = fork,
				extentForward = extent,
				)) - {dnaPoly}

			if nonPolymeraseMolecules:
				print 'Encounter molecules within fork extension range'
				continue

			# Unbind the DNA poly while moving the fork
			self.replicationForksRegions.moleculeLocationIsUnbound(dnaPoly)

			# Move the fork
			newForkPosition = self.replicationForksRegions.forkExtend(
				fork,
				extent- self.dnaPolyForwardFootprint
				)

			# Replace the DNA poly
			self.replicationForksRegions.moleculeLocationIsFork(
				molecule = dnaPoly,
				fork = fork,
				extentForward = self.dnaPolyForwardFootprint,
				extentReverse = self.dnaPolyReverseFootprint
				)

