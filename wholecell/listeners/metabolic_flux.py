#!/usr/bin/env python

"""
MetabolicFlux

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/10/2014
"""

from __future__ import division

import numpy as np
import tables

import wholecell.listeners.listener

# TODO: figure out how to perform non-blocking output

# import matplotlib.pyplot as plt

class MetabolicFlux(wholecell.listeners.listener.Listener):
	""" MetabolicFlux """

	_name = 'MetabolicFlux'

	# Constructor
	def __init__(self, *args, **kwargs):
		# References to other states
		self.states = None

		super(MetabolicFlux, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, kb):
		super(MetabolicFlux, self).initialize(sim, kb)

		self.metabolism = sim.processes['MetabolismFba']

		self.cellCycleLen = kb.cellCycleLen.to('s').magnitude


	# Allocate memory
	def allocate(self):
		super(MetabolicFlux, self).allocate()

		self.fluxes = np.zeros_like(self.metabolism.fluxes)


	def update(self):
		self.fluxes = self.metabolism.fluxes

		# plt.plot(self.fluxes)

	# TODO: implement saving/loading

	# def pytablesCreate(self, h5file, expectedRows):
	# 	# Columns
	# 	d = {
	# 		"time": tables.Int64Col(),
	# 		"cell": tables.Float64Col(),
	# 		"cellDry": tables.Float64Col(),
	# 		"metabolite": tables.Float64Col(),
	# 		"rna": tables.Float64Col(),
	# 		"rrna": tables.Float64Col(),
	# 		"rrna23S": tables.Float64Col(),
	# 		"rrna16S": tables.Float64Col(),
	# 		"rrna5S": tables.Float64Col(),
	# 		"trna": tables.Float64Col(),
	# 		"mrna": tables.Float64Col(),
	# 		"protein": tables.Float64Col(),
	# 		"water": tables.Float64Col(),
	# 		"nucleoid": tables.Float64Col(),
	# 		}

	# 	# Create table
	# 	# TODO: Add compression options (using filters)
	# 	t = h5file.create_table(
	# 		h5file.root,
	# 		self._name,
	# 		d,
	# 		title = self._name,
	# 		filters = tables.Filters(complevel = 9, complib="zlib"),
	# 		expectedrows = expectedRows
	# 		)

	# 	# Store units as metadata
	# 	t.attrs.cell_units = self.massUnits
	# 	t.attrs.cellDry_units = self.massUnits
	# 	t.attrs.metabolite_units = self.massUnits
	# 	t.attrs.rna_units = self.massUnits
	# 	t.attrs.rrna_units = self.massUnits
	# 	t.attrs.rrna23S_units = self.massUnits
	# 	t.attrs.rrna16S_units = self.massUnits
	# 	t.attrs.rrna5S_units = self.massUnits
	# 	t.attrs.trna_units = self.massUnits
	# 	t.attrs.mrna_units = self.massUnits
	# 	t.attrs.protein_units = self.massUnits
	# 	t.attrs.water_units = self.massUnits
	# 	t.attrs.nucleoid_units = self.massUnits


	# def pytablesAppend(self, h5file):
	# 	t = h5file.get_node("/", self._name)
	# 	entry = t.row

	# 	entry["time"] = self.timeStep()
	# 	entry["cell"] = self.cell
	# 	entry["cellDry"] = self.cellDry
	# 	entry["metabolite"] = self.metabolite
	# 	entry["rna"] = self.rna
	# 	entry["rrna"] = self.rrna
	# 	entry["rrna23S"] = self.rrna23S
	# 	entry["rrna16S"] = self.rrna16S
	# 	entry["rrna5S"] = self.rrna5S
	# 	entry["trna"] = self.trna
	# 	entry["mrna"] = self.mrna
	# 	entry["protein"] = self.protein
	# 	entry["water"] = self.water
	# 	entry["nucleoid"] = self.nucleoid

	# 	entry.append()

	# 	t.flush()
