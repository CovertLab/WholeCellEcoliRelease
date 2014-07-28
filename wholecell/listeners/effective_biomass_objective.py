#!/usr/bin/env python

"""
EffectiveBiomassObjective

Records the effective biomass objective as metabolism adapts to the existing 
cell state.

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/24/2014
"""

from __future__ import division

import numpy as np
import tables

import wholecell.listeners.listener

class EffectiveBiomassObjective(wholecell.listeners.listener.Listener):
	""" EffectiveBiomassObjective """

	_name = "EffectiveBiomassObjective"

	# Constructor
	def __init__(self, *args, **kwargs):
		super(EffectiveBiomassObjective, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, kb):
		super(EffectiveBiomassObjective, self).initialize(sim, kb)

		self.effectiveBiomassObjective = None
		self.biomassObjectiveIds = kb.wildtypeBiomass["metaboliteId"]
		self.standardBiomassObjective = kb.wildtypeBiomass["biomassFlux"].to("millimole/DCW_g").magnitude


	# Allocate memory
	def allocate(self):
		super(EffectiveBiomassObjective, self).allocate()

		self.effectiveBiomassObjective = np.zeros_like(self.standardBiomassObjective)
		

	def pytablesCreate(self, h5file, expectedRows):

		# Columns
		dtype = {
			"time": tables.Float64Col(),
			"timeStep": tables.Int64Col(),
			"effectiveBiomassObjective": tables.Float64Col(self.effectiveBiomassObjective.shape)
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

		table.attrs.metaboliteIds = self.biomassObjectiveIds
		table.attrs.standardBiomass = self.standardBiomassObjective


	def pytablesAppend(self, h5file):

		table = h5file.get_node("/", self._name)
		entry = table.row

		entry["time"] = self.time()
		entry["timeStep"] = self.timeStep()
		entry["effectiveBiomassObjective"] = self.effectiveBiomassObjective

		entry.append()

		table.flush()
