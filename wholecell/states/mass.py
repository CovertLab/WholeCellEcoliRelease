#!/usr/bin/env python

"""
Mass

Mass state variable. Represents the total cellular mass.

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/29/2013
"""

from __future__ import division

import numpy as np
import tables

import wholecell.states.state

class Mass(wholecell.states.state.State):
	""" Mass """

	_name = 'Mass'

	# Constructor
	def __init__(self, *args, **kwargs):
		# References to other states
		self.bulkMolecules = None

		# NOTE: molecule weight is converted to femtograms/molecule from
		# grams/mol in BulkMolecules
		self.massUnits = 'fg'

		super(Mass, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, kb):
		super(Mass, self).initialize(sim, kb)

		self.bulkMolecules = sim.states["BulkMolecules"]


	# Allocate memory
	def allocate(self):
		super(Mass, self).allocate()

		self.cell = 0
		self.cellDry = 0
		self.metabolite = 0
		self.rna = 0
		self.protein = 0
		self.nucleoid = 0

		self.growth = 0


	def calculate(self):
		oldMass = self.cell

		self.cell = self.bulkMolecules.mass()

		self.metabolite = self.bulkMolecules.massByType('metabolites')
		self.rna = self.bulkMolecules.massByType('rnas')
		self.protein = self.bulkMolecules.massByType('proteins')
		self.water = self.bulkMolecules.massByType('water')

		self.nucleoid = self.bulkMolecules.massByCompartment('n')

		self.cellDry = self.cell - self.water

		self.growth = self.cell - oldMass


	def pytablesCreate(self, h5file, expectedRows):
		# Columns
		d = {
			"time": tables.Int64Col(),
			"cell": tables.Float64Col(),
			"cellDry": tables.Float64Col(),
			"metabolite": tables.Float64Col(),
			"rna": tables.Float64Col(),
			"protein": tables.Float64Col(),
			"water": tables.Float64Col(),
			"nucleoid": tables.Float64Col(),
			}

		# Create table
		# TODO: Add compression options (using filters)
		t = h5file.create_table(
			h5file.root,
			self._name,
			d,
			title = self._name,
			filters = tables.Filters(complevel = 9, complib="zlib"),
			expectedrows = expectedRows
			)

		# Store units as metadata
		t.attrs.cell_units = self.massUnits
		t.attrs.cellDry_units = self.massUnits
		t.attrs.metabolite_units = self.massUnits
		t.attrs.rna_units = self.massUnits
		t.attrs.protein_units = self.massUnits
		t.attrs.water_units = self.massUnits
		t.attrs.nucleoid_units = self.massUnits


	def pytablesAppend(self, h5file):
		t = h5file.get_node("/", self._name)
		entry = t.row

		entry["time"] = self.timeStep()
		entry["cell"] = self.cell
		entry["cellDry"] = self.cellDry
		entry["metabolite"] = self.metabolite
		entry["rna"] = self.rna
		entry["protein"] = self.protein
		entry["water"] = self.water
		entry["nucleoid"] = self.nucleoid

		entry.append()

		t.flush()
