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
		self.elongationRate = 100

		super(ToyReplication, self).__init__()


	# Construct object graph
	def initialize(self, sim, kb, kb2):
		super(ToyReplication, self).initialize(sim, kb, kb2)

		self.replicationForks = self.chromosomeForksView( # special view for forks
			extentForward = max(self.elongationRate, self.dnaPolyForwardFootprint),
			extentReverse = self.dnaPolyReverseFootprint
			)


	def calculateRequest(self):
		self.replicationForks.requestAll()


	# Calculate temporal evolution
	def evolveState(self):
		# TODO: check whether any forks would collide

		for fork in self.replicationForks.forks():
			# Get DNA polymerases near the fork
			dnaPolymerase = self.replicationForks.moleculeOnFork( # returns the identity of the molecule on the fork, if any
				fork = fork,
				name = 'DNA polymerase'
				)

			# Break if no DNA polymerase
			if dnaPolymerase is None:
				print 'No DNA polymerase on fork'
				
				continue

			# Determine how far we can extend
			extent = self.replicationForks.maximumExtent( # how far we can move without hitting 1) the end of the partitioned space or 2) a fork, up to "extent"
				molecule = dnaPolymerase,
				extentForward = self.elongationRate,
				extentReverse = 0
				)

			# TODO: check/use resources

			# Check for molecules in the way of replication
			nonPolymeraseMolecules = self.replicationForks.moleculesNearFork(
				fork = fork,
				extentForward = extent,
				extentReverse = 0
				) - {dnaPolymerase}

			if nonPolymeraseMolecules:
				# TODO: remove molecules on final position, and move molecules
				# passed onto a random strand in the same position
				print 'Encountered molecules within fork extension range'
				
				continue

			# Unbind the DNA polymerase, for the moment
			self.replicationForks.moleculeLocationIsUnbound(dnaPolymerase)

			# Move the fork
			newForkPosition = self.replicationForks.forkExtend(fork, extent)

			# Place the DNA polymerase back onto the fork
			self.replicationForks.moleculeLocationIsFork( # special moleculeMove routine that places a molecule directly onto a fork
				molecule = dnaPolymerase,
				fork = fork
				)

