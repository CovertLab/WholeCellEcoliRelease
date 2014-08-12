#!/usr/bin/env python

"""
FBAResults

Records dynamics of FBA output.

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/24/2014
"""

from __future__ import division

import numpy as np
import tables

import wholecell.listeners.listener

class FBAResults(wholecell.listeners.listener.Listener):
	""" FBAResults """

	_name = "FBAResults"

	# Constructor
	def __init__(self, *args, **kwargs):
		super(FBAResults, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, kb):
		super(FBAResults, self).initialize(sim, kb)

		self.metabolism = sim.processes["Metabolism"]

		self.objectiveValue = 0.0


	# Allocate memory
	def allocate(self):
		super(FBAResults, self).allocate()

		fba = self.metabolism.fba

		self.reactionIDs = fba.reactionIDs()
		self.externalMoleculeIDs = fba.externalMoleculeIDs()
		self.outputMoleculeIDs = fba.outputMoleculeIDs()

		self.reactionFluxes = np.zeros(len(self.reactionIDs), np.float64)
		self.externalExchangeFluxes = np.zeros(len(self.externalMoleculeIDs), np.float64)
		self.outputFluxes = np.zeros(len(self.outputMoleculeIDs), np.float64)
		self.objectiveComponents = np.zeros_like(self.outputFluxes)


	def pytablesCreate(self, h5file, expectedRows):

		# Columns
		dtype = {
			"time": tables.Float64Col(),
			"timeStep": tables.Int64Col(),
			"reactionFluxes": tables.Float64Col(self.reactionFluxes.shape),
			"externalExchangeFluxes": tables.Float64Col(self.externalExchangeFluxes.shape),
			"outputFluxes": tables.Float64Col(self.outputFluxes.shape),
			"objectiveValue": tables.Float64Col(),
			"objectiveComponents": tables.Float64Col(self.objectiveComponents.shape),
			}

		# Create table
		table = h5file.create_table(
			h5file.root,
			self._name,
			dtype,
			title = self._name,
			filters = tables.Filters(complevel = 9, complib="zlib"),
			expectedrows = expectedRows
			)
	
		groupNames = h5file.create_group(h5file.root,
			'names', 'Reaction and molecule names')

		h5file.create_array(groupNames, 'reactionIDs', [str(s) for s in self.reactionIDs])
		h5file.create_array(groupNames, 'externalMoleculeIDs', [str(s) for s in self.externalMoleculeIDs])
		h5file.create_array(groupNames, 'outputMoleculeIDs', [str(s) for s in self.outputMoleculeIDs])


	def pytablesAppend(self, h5file):

		table = h5file.get_node("/", self._name)
		entry = table.row

		entry["time"] = self.time()
		entry["timeStep"] = self.timeStep()
		entry["reactionFluxes"] = self.reactionFluxes
		entry["externalExchangeFluxes"] = self.externalExchangeFluxes
		entry["outputFluxes"] = self.outputFluxes
		entry["objectiveValue"] = self.objectiveValue
		entry["objectiveComponents"] = self.objectiveComponents

		entry.append()

		table.flush()
