#!/usr/bin/env python
"""
RnapData

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/18/15
"""

from __future__ import division

import numpy as np

import wholecell.listeners.listener

VERBOSE = False


class RnapData(wholecell.listeners.listener.Listener):
	""" RnapData """

	_name = 'RnapData'

	# Constructor
	def __init__(self, *args, **kwargs):
		super(RnapData, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, sim_data):
		super(RnapData, self).initialize(sim, sim_data)

		self.nRnaSpecies = sim_data.process.transcription.rnaData['id'].size
		self.uniqueMolecules = sim.internal_states['UniqueMolecules']


	# Allocate memory
	def allocate(self):
		super(RnapData, self).allocate()

		# Positions of active RNAPs on the chromosome
		# The size of this array must be larger than the maximum possible
		# counts of active RNAPs at any timestep of the simulation.
		self.active_rnap_coordinates = np.full(10000, np.nan, np.float64)

		# Attributes broadcast by the PolypeptideElongation process
		self.actualElongations = 0
		self.didTerminate = 0
		self.didInitialize = 0
		self.terminationLoss = 0
		self.rnaInitEvent = np.zeros(self.nRnaSpecies, np.int64)

		# Collisions with replisomes
		self.n_aborted_initiations = 0
		self.n_collisions = 0


	def update(self):
		self.active_rnap_coordinates[:] = np.nan

		active_rnaps = self.uniqueMolecules.container.objectsInCollection(
			'activeRnaPoly')

		# Read coordinates of all active RNAPs
		if len(active_rnaps) > 0:
			active_rnap_coordinates = active_rnaps.attr("coordinates")
			self.active_rnap_coordinates[
			:active_rnap_coordinates.size] = active_rnap_coordinates


	def tableCreate(self, tableWriter):
		pass


	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			active_rnap_coordinates=self.active_rnap_coordinates,
			actualElongations = self.actualElongations,
			didTerminate = self.didTerminate,
			didInitialize = self.didInitialize,
			terminationLoss = self.terminationLoss,
			rnaInitEvent = self.rnaInitEvent,
			n_aborted_initiations=self.n_aborted_initiations,
			n_collisions=self.n_collisions,
			)
