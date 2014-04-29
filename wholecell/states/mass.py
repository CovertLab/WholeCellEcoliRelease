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

STATES_WITH_MASS = ['BulkMolecules', 'UniqueMolecules']

class Mass(wholecell.states.state.State):
	""" Mass """

	_name = 'Mass'

	# Constructor
	def __init__(self, *args, **kwargs):
		# References to other states
		self.states = None

		# NOTE: molecule weight is converted to femtograms/molecule from
		# grams/mol in BulkMolecules
		self.massUnits = 'fg'

		super(Mass, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, kb):
		super(Mass, self).initialize(sim, kb)

		self.bulkMolecules = sim.states['BulkMolecules']

		self.states = sim.states

		self.cellCycleLen = kb.cellCycleLen.to('s').magnitude


	# Allocate memory
	def allocate(self):
		super(Mass, self).allocate()

		self.setInitial = False

		self.cellDryInitial = None
		self.proteinInitial = None
		self.rnaInitial = None

		self._resetMasses()

		self.cellDry = 0

		self.proteinFraction = 0
		self.rnaFraction = 0
		
		self.cellDryFoldChange = 0
		self.proteinFoldChange = 0
		self.rnaFoldChange = 0

		self.expectedFoldChange = 0

		self.growth = 0


	def _resetMasses(self):
		self.cell = 0

		self.metabolite = 0
		self.rna = 0
		self.rrna = 0
		self.protein = 0
		self.nucleoid = 0

		self.water = 0


	def calculate(self):
		oldMass = self.cellDry

		self._resetMasses()

		# TODO: rework mass calculations as State methods
		# TODO: change states without mass (such as Mass itself) into "listener" classes
		for stateName, state in self.states.viewitems():
			if stateName not in STATES_WITH_MASS: continue

			self.cell += state.mass()

			self.metabolite += state.massByType('metabolites')
			self.rna += state.massByType('rnas')
			self.protein += state.massByType('proteins')
			self.nucleoid += state.massByCompartment('n')
			
			self.water += state.massByType('water')
		
		self.rrna = self.bulkMolecules.massByType('rrnas')
		self.rrna23S = self.bulkMolecules.massByType('rrna23Ss')
		self.rrna16S = self.bulkMolecules.massByType('rrna16Ss')
		self.rrna5S = self.bulkMolecules.massByType('rrna5Ss')
		self.trna = self.bulkMolecules.massByType('trnas')
		self.mrna = self.bulkMolecules.massByType('mrnas')

		self.cellDry = self.cell - self.water

		self.growth = self.cellDry - oldMass

		self.proteinFraction = self.protein / self.cellDry
		self.rnaFraction = self.rna / self.cellDry

		if not self.setInitial:
			self.setInitial = True

			self.cellDryInitial = self.cellDry
			self.proteinInitial = self.protein
			self.rnaInitial = self.rna
			self.rrnaInitial = self.rrna
			self.rrna23SInitial = self.rrna23S
			self.rrna16SInitial = self.rrna16S
			self.rrna5SInitial = self.rrna5S
			self.trnaInitial = self.trna
			self.mrnaInitial = self.mrna

		self.cellDryFoldChange = self.cellDry / self.cellDryInitial
		self.proteinFoldChange = self.protein / self.proteinInitial
		self.rnaFoldChange = self.rna / self.rnaInitial
		self.rrnaFoldChange = self.rrna / self.rrnaInitial
		self.rrna23SFoldChange = self.rrna23S / self.rrna23SInitial
		self.rrna16SFoldChange = self.rrna16S / self.rrna16SInitial
		self.rrna5SFoldChange = self.rrna5S / self.rrna5SInitial
		self.trnaFoldChange = self.trna / self.trnaInitial
		self.mrnaFoldChange = self.mrna / self.mrnaInitial

		self.expectedFoldChange = np.exp(np.log(2) * self.time() / self.cellCycleLen)

		# print (
		# 	self.rnaFoldChange, self.rrnaFoldChange, self.trnaFoldChange, self.mrnaFoldChange,
		# 	self.rrna23SFoldChange, self.rrna16SFoldChange, self.rrna5SFoldChange
		# 	)


	def pytablesCreate(self, h5file, expectedRows):
		# Columns
		d = {
			"time": tables.Int64Col(),
			"cell": tables.Float64Col(),
			"cellDry": tables.Float64Col(),
			"metabolite": tables.Float64Col(),
			"rna": tables.Float64Col(),
			"rrna": tables.Float64Col(),
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
		t.attrs.rrna_units = self.massUnits
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
		entry["rrna"] = self.rrna
		entry["protein"] = self.protein
		entry["water"] = self.water
		entry["nucleoid"] = self.nucleoid

		entry.append()

		t.flush()
