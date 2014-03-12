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

		self.dnaPolyForwardFootprint = 50
		self.dnaPolyReverseFootprint = 50
		self.elongationRate = 100

		super(ToyReplication, self).__init__()


	# Construct object graph
	def initialize(self, sim, kb):
		super(ToyReplication, self).initialize(sim, kb)

		self.chromosome = self.chromosomeView() # a view which can operate on the partitioned state but does not run requests/queries

		self.replicationForks = self.replicationForkView( # a view which can request/query forks and return locations of forks on the partitioned state
			extentForward = max(self.elongationRate, self.dnaPolyForwardFootprint),
			extentReverse = self.dnaPolyReverseFootprint
			)


	def calculateRequest(self):
		self.replicationForks.requestAll()


	# Calculate temporal evolution
	def evolveState(self):
		# TODO: check whether any forks would collide

		for forkStrand, forkPosition, forkDirection in self.replicationForks.forks():
			# Get DNA polymerases near the fork
			dnaPolymerase = self.chromosome.moleculeOnFork( # returns the identity of the molecule on the fork, if any
				name = 'DNA polymerase'
				forkStrand = forkStrand,
				forkPosition = forkPosition
				)

			# Break if no DNA polymerase
			if dnaPolymerase is None:
				print 'No DNA polymerase on fork'
				break

			# Determine how far we can extend
			extent = self.chromosome.maximumExtent( # how far we can move without hitting 1) the end of the partitioned space or 2) a fork, up to "extent"
				strand = forkStrand,
				position = forkPosition,
				direction = forkDirection,
				extent = self.elongationRate
				)

			# TODO: check/use resources

			nonPolymeraseMolecules = self.chromosome.molecules( # returns molecules that satisfy the arguments
				strand = forkStrand,
				position = forkPosition,
				direction = forkDirection,
				extent = extent
				) - {dnaPolymerase}

			if nonPolymeraseMolecules:
				# TODO: remove molecules on final position, and move molecules
				# passed onto a random strand in the same position
				print 'Encountered molecules within fork extension range'
				break

			newForkPosition = self.chromosome.extendFork(forkStrand, # extends an indicated fork (direction is inferred)
				forkPosition, extent)

			self.chromosome.moveMoleculeToFork( # special moveMolecule routine that places a molecule directly onto a fork
				dnaPolymerase,
				forkStrand = forkStrand,
				forkPosition = newForkPosition,
				extentForward = self.dnaPolyForwardFootprint,
				extentReverse = self.dnaPolyReverseFootprint
				)

