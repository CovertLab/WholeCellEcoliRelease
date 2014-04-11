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

	compartments = [
		{"id": "c", "name": "Cytosol"},
		{"id": "e", "name": "Extracellular space"},
		{"id": "i", "name": "Inner membrane"},
		{"id": "j", "name": "Projection"},
		{"id": "l", "name": "Pilus"},
		{"id": "m", "name": "Membrane"},
		{"id": "n", "name": "Nucleoid"},
		{"id": "o", "name": "Outer membrane"},
		{"id": "p", "name": "Periplasm"},
		{"id": "w", "name": "Cell wall"}
		]
	
	cIdx = {c['id']:i for i, c in enumerate(compartments)}

	# Constructor
	def __init__(self, *args, **kwargs):
		self.meta = {
			"id": "Mass",
			"name": "Mass",
			"dynamics": ["total", "cell", "cellDry", "metabolite", "rna", "protein", 'growth'],
			"units": {
				"total": "fg",
				"cellDry": "fg",
				"metabolite": "fg",
				"rna": "fg",
				"protein": "fg",
				'growth': "fg/s"
				}
		}

		# References to other states
		self.bulkMolecules = None
		self.time = None

		# Mass
		self.cell = None
		self.cellDry = None
		self.metabolite = None
		self.rna = None
		self.protein = None

		self.growth = None

		super(Mass, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, kb):
		super(Mass, self).initialize(sim, kb)

		self.bulkMolecules = sim.states["BulkMolecules"]
		self.time = sim.states["Time"]

		self.nAvogadro = kb.nAvogadro.to('1/mole').magnitude


	# Allocate memory
	def allocate(self):
		super(Mass, self).allocate()

		# TODO: reimplement compartment-specific records, if desired

		self.cell = np.zeros(1)
		self.cellDry = np.zeros(1)
		self.metabolite = np.zeros(1)
		self.rna = np.zeros(1)
		self.protein = np.zeros(1)

		self.growth = np.zeros(1)


	def calculate(self):
		oldMass = self.cell

		# Total
		self.cell = self.bulkMolecules.mass() / self.nAvogadro * 1e15

		# Cell
		self.metabolite = self.bulkMolecules.mass('metabolites') / self.nAvogadro * 1e15
		self.rna        = self.bulkMolecules.mass('rnas')        / self.nAvogadro * 1e15
		self.protein    = self.bulkMolecules.mass('proteins')    / self.nAvogadro * 1e15

		self.cellDry = self.cell - self.bulkMolecules.mass('water') / self.nAvogadro * 1e15

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
			}

		# Create table
		# TODO: Add compression options (using filters)
		t = h5file.create_table(
			h5file.root,
			self.meta["id"],
			d,
			title = self.meta["name"],
			filters = tables.Filters(complevel = 9, complib="zlib"),
			expectedrows = expectedRows
			)

		# Store units as metadata
		t.attrs.cell_units = self.meta["units"]["cell"]
		t.attrs.cellDry_units = self.meta["units"]["cellDry"]
		t.attrs.metabolite_units = self.meta["units"]["metabolite"]
		t.attrs.rna_units = self.meta["units"]["rna"]
		t.attrs.protein_units = self.meta["units"]["protein"]


	def pytablesAppend(self, h5file):
		simTime = self.time.value
		t = h5file.get_node("/", self.meta["id"])
		entry = t.row

		entry["time"] = simTime
		entry["cell"] = self.cell
		entry["cellDry"] = self.cellDry
		entry["metabolite"] = self.metabolite
		entry["rna"] = self.rna
		entry["protein"] = self.protein

		entry.append()

		t.flush()
