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

	meta = {
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

	compartments = [
		{"id": "c", "name": "Cytosol"},
		{"id": "e", "name": "Extracellular space"},
		{"id": "m", "name": "Membrane"}
	]
	cIdx = {"c": 1, "e": 2, "m": 3}

	# Constructor
	def __init__(self, *args, **kwargs):
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
		super(Time, self).allocate()

		self.total = numpy.zeros((1, len(self.compartments)))
		self.cell = numpy.zeros((1, len(self.compartments)))
		self.cellDry = numpy.zeros((1, len(self.compartments)))
		self.metabolite = numpy.zeros((1, len(self.compartments)))
		self.rna = numpy.zeros((1, len(self.compartments)))
		self.protein = numpy.zeros((1, len(self.compartments)))

	# Calculate (and cache) any dependent properties
	def calculate(self):
		from wholecell.sim.util.Constants import Constants

		mc = self.moleculeCounts

		# Total
		self.total = ( numpy.dot(mc.mws, mc.counts) ) / Constants.nAvogadro * 1e15

		# Cell
		self.metabolite = ( numpy.dot(mc.mws[mc.types == mc.typeVals.metabolite], mc.counts[mc.types == mc.typeVals.metabolite]) ) / Constants.nAvogadro * 1e15
		self.rna        = ( numpy.dot(mc.mws[mc.types == mc.typeVals.rna       ], mc.counts[mc.types == mc.typeVals.rna       ]) ) / Constants.nAvogadro * 1e15
		self.protein    = ( numpy.dot(mc.mws[mc.types == mc.typeVals.protein   ], mc.counts[mc.types == mc.typeVals.protein   ]) ) / Constants.nAvogadro * 1e15

		cIdxs = numpy.array([self.cIdx.c, self.cIdx.m])

		self.cell[:] = 0
		this.cell[cIdxs] = self.metabolite[cIdxs] + self.rna[cIdxs] + self.protein[cIdxs]

		self.cellDry[:] = 0
		self.cellDry[cIdxs] = self.cell[cIdxs] - ( numpy.dot(mc.mws[mc.idx.h2o], mc.counts[mc.idx.h2o, cIdxs]) ) / Constants.nAvogadro * 1e15