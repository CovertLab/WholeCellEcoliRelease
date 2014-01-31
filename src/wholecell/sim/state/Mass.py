#!/usr/bin/env python

"""
Mass

Mass state variable. Represents the total cellular mass.

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/29/2013
"""

import numpy

import wholecell.sim.state.State

class Mass(wholecell.sim.state.State.State):
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
	cIdx = {"c": 0, "e": 1, "i": 2, "j": 3, "l": 4, "m": 5, "n": 6, "o": 7, "p": 8, "w": 9}

	# Constructor
	def __init__(self, *args, **kwargs):
		self.meta = {
			"id": "Mass",
			"name": "Mass",
			"dynamics": ["total", "cell", "cellDry", "metabolite", "rna", "protein", 'growth'],
			"units": {
				"total": "fg",
				"cell": "fg",
				"cellDry": "fg",
				"metabolite": "fg",
				"rna": "fg",
				"protein": "fg",
				'growth': "fg/s"
				}
		}

		# References to other states
		self.moleculeCounts = None
		self.time = None

		# Mass
		self.total = None
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

		self.moleculeCounts = sim.getState("MoleculeCounts")
		self.time = sim.getState("Time")

	# Allocate memory
	def allocate(self):
		super(Mass, self).allocate()

		self.total = numpy.zeros(len(self.compartments))
		self.cell = numpy.zeros(len(self.compartments))
		self.cellDry = numpy.zeros(len(self.compartments))
		self.metabolite = numpy.zeros(len(self.compartments))
		self.rna = numpy.zeros(len(self.compartments))
		self.protein = numpy.zeros(len(self.compartments))

		self.growth = numpy.zeros(1)


	def calculate(self):
		from wholecell.util.Constants import Constants

		mc = self.moleculeCounts

		# Total
		self.total = mc.massAll() / Constants.nAvogadro * 1e15

		# Cell
		self.metabolite = mc.massAll('metabolites') / Constants.nAvogadro * 1e15
		self.rna        = mc.massAll('rnas')        / Constants.nAvogadro * 1e15
		self.protein    = mc.massAll('proteins')    / Constants.nAvogadro * 1e15

		cIdxs = numpy.array([
							self.cIdx["c"], self.cIdx["i"], self.cIdx["j"], self.cIdx["l"], self.cIdx["m"],
							self.cIdx["n"], self.cIdx["o"], self.cIdx["p"], self.cIdx["w"]])

		oldMass = self.cell.sum()

		self.cell[:] = 0
		self.cell[cIdxs] = self.metabolite[cIdxs] + self.rna[cIdxs] + self.protein[cIdxs]

		self.cellDry[:] = 0
		self.cellDry[cIdxs] = self.cell[cIdxs] - mc.massAll('water')[cIdxs] / Constants.nAvogadro * 1e15

		self.growth = self.cell.sum() - oldMass


	def pytablesCreate(self, h5file):
		import tables

		# Columns
		d = {
			"time": tables.Int64Col(),
			"compartment": tables.StringCol(max([len(x) for x in self.cIdx.keys()])),
			"total": tables.Float64Col(),
			"cell": tables.Float64Col(),
			"cellDry": tables.Float64Col(),
			"metabolite": tables.Float64Col(),
			"rna": tables.Float64Col(),
			"protein": tables.Float64Col(),
			}

		# Create table
		# TODO: Add compression options (using filters)
		t = h5file.create_table(h5file.root, self.meta["id"], d, title = self.meta["name"], filters = tables.Filters(complevel = 9, complib="zlib"))

		# Store units as metadata
		t.attrs.total_units = self.meta["units"]["total"]
		t.attrs.cell_units = self.meta["units"]["cell"]
		t.attrs.cellDry_units = self.meta["units"]["cellDry"]
		t.attrs.metabolite_units = self.meta["units"]["metabolite"]
		t.attrs.rna_units = self.meta["units"]["rna"]
		t.attrs.protein_units = self.meta["units"]["protein"]

	def pytablesAppend(self, h5file):
		import tables

		simTime = self.time.value
		t = h5file.get_node("/", self.meta["id"])
		entry = t.row

		for i in xrange(len(self.cIdx.keys())):
				entry["time"] = simTime
				entry["compartment"] = [key for key,val in self.cIdx.iteritems() if val == i][0]
				entry["total"] = self.total[i]
				entry["cell"] = self.cell[i]
				entry["cellDry"] = self.cellDry[i]
				entry["metabolite"] = self.metabolite[i]
				entry["rna"] = self.rna[i] # Ha! RNAi!
				entry["protein"] = self.protein[i]
				entry.append()

		t.flush()
