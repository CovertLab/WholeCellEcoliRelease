#!/usr/bin/env python

"""
Metabolism

Metabolism sub-model. Encodes molecular simulation of microbial metabolism using flux-balance analysis.

TODO:
- move over to flexFBA
- implement metabolite pools
- enzyme-limited reactions (& fit enzyme expression)
- option to call a reduced form of metabolism (assume optimal)

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/2/2013
"""

from __future__ import division

import numpy as np

import wholecell.processes.process

from wholecell.utils.constants import REQUEST_PRIORITY_METABOLISM

class Metabolism(wholecell.processes.process.Process):
	""" Metabolism """

	_name = "Metabolism"

	# Constructor
	def __init__(self):
		super(Metabolism, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(Metabolism, self).initialize(sim, kb)

		self.fba = FluxBalanceAnalysis('arguments...')

		self.fba.externalMoleculeCounts('max # of molecules that can be used in a time step, defined by media/diffusion')

		self.molecules = self.bulkMoleculesView('internal molecule names')

		self.enzymes = self.bulkMoleculesView('enzyme names')

		self.bulkMoleculeRequestPriorityIs(REQUEST_PRIORITY_METABOLISM)


	def calculateRequest(self):
		self.molecules.requestAll()


	# Calculate temporal evolution
	def evolveState(self):
		# Setup
		self.fba.internalMoleculeCountsIs(self.molecules.counts())
		self.fba.enzymeCountsIs(self.enzymes.counts())

		# Run
		self.fba.run()

		# Finalize
		self.molecules.countsIs(self.fba.internalMoleculeCounts())

		# TODO: record solution metadata, probably in a listener


class FluxBalanceAnalysis(object):
	""" FluxBalanceAnalysis

	Solver for various FBA implementations.

	
	Required arguments:

	- reactionData, a dict of strings:dicts (reactionID:reaction information dict)
		Each value in the dict is another dict with items
		- enzymeID : string or None
		- stoichiometry : dict of molecule ID to stoichiometry pairs
		- isReversible : bool

	- externalExchangedMolecules, an iterable of strings (moleculeIDs)
		Every provided ID will be set up with an exchange flux.

	- fluxObjective, a dict of strings:floats (metaboliteId:objective value)
		The meaning and usage of the objective will vary depending on the 
		formulation of FBA desired.

	
	Optional arguments (set to None for default behavior):

	- objectiveType, a string
		"standard": standard FBA objective (default)
		"flexible": flexFBA
		"pools": similar to FBA; optimizes towards desired pool concentrations

	- objectiveParameters, a dict
		Keys are specific to the FBA implementation.

	- internalExchangedMolecules, an iterable of strings (moleculeIDs)
		Every provided ID will be set up with an exchange flux.
	
	- enzymeRates, a dict of strings:floats (enzymeID:catalytic rate constant * dt)
		Used to set up the pseudo metabolites and boundary constraints needed
		for limiting reactions by enzyme counts.

	- moleculeMasses, a dict of floats (moleculeID:mass)
		Used in computing the net mass into the system.  Only needed and used 
		for moleculeIDs in externalExchangedMolecules.

	
	Caveats:
	
	There is no strict type checking, despite what the above may imply.

	During initialization, an exception will be raised if a reference is made 
	to an unutilized metabolite/enzyme/etc as described by the reaction 
	network.

	"""

	def __init__(self, reactionData, externalExchangedMolecules, fluxObjective,
			objectiveType = None, objectiveParameters = None,
			internalExchangedMolecules = None, enzymeRates = None,
			moleculeMasses = None):

		# Set up lists of node (row) and edge (column) names

		# Set up IJV lists for sparse array

		# Parses reactions into molecule and reaction names, and add stoich
		# to IJV lists

		# Add exchangedMoleculeIDs to IJV lists

		# Set up the objective (implementation varies)

		# Set up enzyme pseudometabolites and constraints

		# Set up mass accumulation column

		# Create cvxopt abstractions

		pass


	# Constraint setup

	def externalMoleculeCountsIs(self, counts):
		pass


	def internalMoleculeCountsIs(self, counts):
		pass


	def enzymeCountsIs(self, counts):
		pass


	# Evaluation

	def run(self):
		pass


	# Output

	def internalMoleculeCounts(self):
		pass

