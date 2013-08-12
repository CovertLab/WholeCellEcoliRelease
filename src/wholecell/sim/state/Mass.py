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
			"dynamics": ["total", "cell", "cellDry", "metabolite", "rna", "protein"],
			"units": {
				"total": "fg",
				"cell": "fg",
				"cellDry": "fg",
				"metabolite": "fg",
				"rna": "fg",
				"protein": "fg"
				}
		}

		# References to other states
		self.moleculeCounts = None

		# Mass
		self.total = None
		self.cell = None
		self.cellDry = None
		self.metabolite = None
		self.rna = None
		self.protein = None

		super(Mass, self).__init__(*args, **kwargs)

	# Construct object graph
	def initialize(self, sim, kb):
		super(Mass, self).initialize(sim, kb)

		self.moleculeCounts = sim.getState("MoleculeCounts")

	# Allocate memory
	def allocate(self):
		super(Mass, self).allocate()

		self.total = numpy.zeros(len(self.compartments))
		self.cell = numpy.zeros(len(self.compartments))
		self.cellDry = numpy.zeros(len(self.compartments))
		self.metabolite = numpy.zeros(len(self.compartments))
		self.rna = numpy.zeros(len(self.compartments))
		self.protein = numpy.zeros(len(self.compartments))

	# Calculate (and cache) any dependent properties
	def calculate(self):
		from wholecell.util.Constants import Constants

		mc = self.moleculeCounts

		# Total
		self.total = ( numpy.dot(mc.mws, mc.counts) ) / Constants.nAvogadro * 1e15

		# Cell
		self.metabolite = ( numpy.dot(mc.mws[mc.types == mc.typeVals["metabolite"]], mc.counts[mc.types == mc.typeVals["metabolite"]]) ) / Constants.nAvogadro * 1e15
		self.rna        = ( numpy.dot(mc.mws[mc.types == mc.typeVals["rna"]       ], mc.counts[mc.types == mc.typeVals["rna"]       ]) ) / Constants.nAvogadro * 1e15
		self.protein    = ( numpy.dot(mc.mws[mc.types == mc.typeVals["protein"]   ], mc.counts[mc.types == mc.typeVals["protein"]   ]) ) / Constants.nAvogadro * 1e15

		cIdxs = numpy.array([
							self.cIdx["c"], self.cIdx["i"], self.cIdx["j"], self.cIdx["l"], self.cIdx["m"],
							self.cIdx["n"], self.cIdx["o"], self.cIdx["p"], self.cIdx["w"]])

		self.cell[:] = 0
		self.cell[cIdxs] = self.metabolite[cIdxs] + self.rna[cIdxs] + self.protein[cIdxs]

		self.cellDry[:] = 0
		self.cellDry[cIdxs] = self.cell[cIdxs] - ( mc.mws[mc.idx["h2o"]] * mc.counts[mc.idx["h2o"], cIdxs] ) / Constants.nAvogadro * 1e15

	def pytablesCreate(self, h5file, sim):
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

	def pytablesAppend(self, h5file, sim):
		import tables

		simTime = sim.getState("Time").value
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