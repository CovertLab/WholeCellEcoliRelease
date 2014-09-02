#!/usr/bin/env python

"""
Mass

Mass listener. Represents the total cellular mass.

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/29/2013
"""

# TODO: generalize this logic for use with a generic simulation

from __future__ import division

import numpy as np
import tables

import wholecell.listeners.listener
from wholecell.utils import units

class Mass(wholecell.listeners.listener.Listener):
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

		self.states = sim.states

		self.processNames = list(sim.processes.keys()) + ["Unallocated"]

		self.cellCycleLen = kb.cellCycleLen.asNumber(units.s)

		self.rnaIndexes = np.array([
			kb.submassNameToIndex[name]
			for name in ["23srRNA", "16srRNA", "5srRNA", "tRNA", "mRNA", "miscRNA", "RNA"]
			])

		self.rRnaIndexes = np.array([
			kb.submassNameToIndex[name]
			for name in ["23srRNA", "16srRNA", "5srRNA"]
			])

		self.tRnaIndex = kb.submassNameToIndex["tRNA"]

		self.mRnaIndex = kb.submassNameToIndex["mRNA"]

		self.dnaIndex = kb.submassNameToIndex["DNA"]

		self.proteinIndex = kb.submassNameToIndex["protein"]

		self.waterIndex = kb.submassNameToIndex["water"]

		# Set initial values

		self.setInitial = False

		self.dryMass = 0
		# TODO: set initial masses based on some calculations of the expected
		# mother cell (divided by two) in the last time step

		# Register logged quantities

		self.registerLoggedQuantity(
			"Cell mass\n(fg)",
			"cellMass",
			".2f"
			)

		self.registerLoggedQuantity(
			"Dry mass\n(fg)",
			"dryMass",
			".2f"
			)

		self.registerLoggedQuantity(
			"Dry mass\nfold change",
			"dryMassFoldChange",
			".3f"
			)

		self.registerLoggedQuantity(
			"Expected\nfold change",
			"expectedMassFoldChange",
			".3f"
			)

		self.registerLoggedQuantity(
			"Growth\n(fg/s)",
			"growth",
			".4f"
			)

		self.registerLoggedQuantity(
			"Protein\nfraction",
			"proteinMassFraction",
			".3f"
			)

		self.registerLoggedQuantity(
			"Protein\nfold change",
			"proteinMassFoldChange",
			".3f"
			)

		self.registerLoggedQuantity(
			"RNA\nfraction",
			"rnaMassFraction",
			".3f"
			)

		self.registerLoggedQuantity(
			"RNA\nfold change",
			"rnaMassFoldChange",
			".3f"
			)


	def update(self):
		oldDryMass = self.dryMass

		masses = sum(state.mass() for state in self.states.itervalues())

		preEvolveMasses = masses[0, ...]
		postEvolveMasses = masses[1, ...]

		self.cellMass = postEvolveMasses.sum() # sum over all dimensions
		submasses = postEvolveMasses.sum(axis = 0) # sum over the processes

		self.waterMass = submasses[self.waterIndex]
		self.dryMass = self.cellMass - self.waterMass
		self.rnaMass = submasses[self.rnaIndexes].sum()
		self.rRnaMass = submasses[self.rRnaIndexes].sum()
		self.tRnaMass = submasses[self.tRnaIndex]
		self.mRnaMass = submasses[self.mRnaIndex]
		self.dnaMass = submasses[self.dnaIndex]
		self.proteinMass = submasses[self.proteinIndex]

		processInitialMass = preEvolveMasses.sum(axis = 1)
		processFinalMass = postEvolveMasses.sum(axis = 1)

		self.processMassDifferences = processFinalMass - processInitialMass

		if self.timeStep() > 0:
			self.growth = self.dryMass - oldDryMass

		else:
			self.growth = 0
			# TODO: solve for an expected initial growth rate

		self.proteinMassFraction = self.proteinMass / self.dryMass
		self.rnaMassFraction = self.rnaMass / self.dryMass

		if not self.setInitial:
			self.setInitial = True

			self.dryMassInitial = self.dryMass
			self.proteinMassInitial = self.proteinMass
			self.rnaMassInitial = self.rnaMass

		self.dryMassFoldChange = self.dryMass / self.dryMassInitial
		self.proteinMassFoldChange = self.proteinMass / self.proteinMassInitial
		self.rnaMassFoldChange = self.rnaMass / self.rnaMassInitial

		self.expectedMassFoldChange = np.exp(np.log(2) * self.time() / self.cellCycleLen)


	def pytablesCreate(self, h5file, expectedRows):
		# Columns
		d = {
			"time": tables.Float64Col(),
			"timeStep": tables.Int64Col(),
			"cellMass": tables.Float64Col(),
			"growth" : tables.Float64Col(),
			"dryMass": tables.Float64Col(),
			"rnaMass": tables.Float64Col(),
			"rRnaMass": tables.Float64Col(),
			"tRnaMass": tables.Float64Col(),
			"mRnaMass": tables.Float64Col(),
			"dnaMass": tables.Float64Col(),
			"proteinMass": tables.Float64Col(),
			"waterMass": tables.Float64Col(),
			"processMassDifferences": tables.Float64Col(len(self.processNames)),
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
		t.attrs.growth_units = self.massUnits
		t.attrs.rna_units = self.massUnits
		t.attrs.protein_units = self.massUnits
		t.attrs.water_units = self.massUnits
		t.attrs.nucleoid_units = self.massUnits

		t.attrs.processNames = self.processNames


	def pytablesAppend(self, h5file):
		t = h5file.get_node("/", self._name)
		entry = t.row

		entry["time"] = self.time()
		entry["timeStep"] = self.timeStep()
		entry["cellMass"] = self.cellMass
		entry["growth"] = self.growth
		entry["dryMass"] = self.dryMass
		entry["rnaMass"] = self.rnaMass
		entry["rRnaMass"] = self.rRnaMass
		entry["tRnaMass"] = self.tRnaMass
		entry["mRnaMass"] = self.mRnaMass
		entry["dnaMass"] = self.dnaMass
		entry["proteinMass"] = self.proteinMass
		entry["waterMass"] = self.waterMass
		entry["processMassDifferences"] = self.processMassDifferences

		entry.append()

		t.flush()
