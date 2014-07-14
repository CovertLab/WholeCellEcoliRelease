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

		self.fba.externalMoleculeCountsIs('max # of molecules that can be used in a time step, defined by media/diffusion')

		self.molecules = self.bulkMoleculesView(self.fba.outputMoleculeIDs())

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
		self.molecules.countsIs(self.fba.outputMoleculeCounts())

		# TODO: record solution metadata, probably in a listener


# TODO: move below to a new file

# import numpy as np
import cvxopt

class FBAException(Exception):
	pass


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

	- objective, a dict of strings:floats (moleculeID:objective value)
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

	# Initialization

	def __init__(self, reactionData, externalExchangedMolecules, objective,
			objectiveType = None, objectiveParameters = None,
			internalExchangedMolecules = None, enzymeRates = None,
			moleculeMasses = None):

		# Set up attributes

		self._nodeNames = []
		self._edgeNames = []

		# Set up running values for initialization

		## Used for creating the sparse matrix
		rowIndexes = []
		colIndexes = []
		values = []

		## Map reactions to some limitations (enzymatic, thermodynamic)
		reactionIndexToEnzyme = {}
		reactionIndexes = []
		reactionIsReversible = []

		## Output calculations
		outputMoleculeIndexes = []
		outputReactionIndexes = []

		# Parses reaction data

		for reactionID, data in reactionData.viewitems():
			colIndex = self._edgeAdd(reactionID)

			for moleculeID, stoichCoeff in data["stoichiometry"].viewitems():
				rowIndex = self._nodeIndex(moleculeID, True)

				rowIndexes.append(rowIndex)
				colIndexes.append(colIndex)
				values.append(stoichCoeff)

			reactionIndexes.append(colIndex)

			if data["enzymeID"] is not None:
				reactionIndexToEnzyme[colIndex] = data["enzyme"]

			reactionIsReversible.append(data["isReversible"])

		self._reactionIndexes = np.array(reactionIndexes)
		self._reactionIsReversible = np.array(reactionIsReversible)

		# Add external exchange reactions

		externalExchangeIndexes = []

		for moleculeID in externalExchangedMolecules:
			exchangeID = "{} external exchange".format(moleculeID)

			colIndex = self._edgeAdd(exchangeID)
			rowIndex = self._nodeIndex(moleculeID)

			# NOTE: The convention, if I am not mistaken, is to define 
			# exchanges as pointing out of the system, even though it is common
			# to think of exchanges as inputs.  Regardless this choice should 
			# have no impact outside of this class.

			rowIndexes.append(rowIndex)
			colIndexes.append(colIndex)
			values.append(-1)

			externalExchangeIndexes.append(colIndex)

		self._externalExchangeIndexes = np.array(externalExchangeIndexes)

		# Set up the objective

		## First, set up the molecule count -> objective equivalents flux

		objectiveEquivalentIndexes = []

		for moleculeID, coeff in objective.viewitems():
			molecule_rowIndex = self._nodeIndex(moleculeID)

			pseudoFluxID = "molecules of {} to fractional objective equivalents".format(moleculeID)
			colIndex = self._edgeAdd(pseudoFluxID)

			objectiveEquivID = "fractional objective equivalent for {}".format(moleculeID)
			objectiveEquiv_rowIndex = self._nodeAdd(objectiveEquivID)

			rowIndexes.append(molecule_rowIndex)
			colIndexes.append(colIndex)
			values.append(-coeff)

			rowIndexes.append(objectiveEquiv_rowIndex)
			colIndexes.append(colIndex)
			values.append(+1)

			outputMoleculeIndexes.append(molecule_rowIndex)
			objectiveEquivalentIndexes.append(objectiveEquiv_rowIndex)

			outputReactionIndexes.append(colIndex)

		## Next, set up the actual objective function (implementation varies)
		objIndexes = []
		objValues = []

		if objectiveType is None or objectiveType == "standard":
			objectiveID = "Standard FBA objective reaction"

			colIndex = self._edgeAdd(objectiveID)

			nObjectiveEquivalents = len(objectiveEquivalentIndexes)

			rowIndexes.extend(objectiveEquivalentIndexes)
			colIndexes.extend([colIndex]*nObjectiveEquivalents)
			values.extend([-1]*nObjectiveEquivalents)

			objIndexes.append(colIndex)
			objValues.append(+1)

		elif objectiveType == "flexible":
			raise NotImplementedError()

		elif objectiveType == "pools":
			raise NotImplementedError()

		else:
			raise FBAException("Unrecognized objectiveType: {}".format(objectiveType))

		# Add internal exchange reactions

		if internalExchangedMolecules is not None:
			raise NotImplementedError()

		# Set up enzyme pseudometabolites and constraints

		if enzymeRates is not None:
			raise NotImplementedError()

		# Set up mass accumulation column

		if moleculeMasses is not None:
			raise NotImplementedError()

		# Finalize some running values

		self._nEdges = len(self._edgeNames)
		self._nNodes = len(self._nodeNames)

		self._outputMoleculeIDs = tuple([self._nodeNames[index] for index in outputMoleculeIndexes])

		self._outputReactionIndexes = np.array(outputReactionIndexes)

		# Create cvxopt abstractions

		# Optimization problem:
		#
		# \max_v f^T x
		# 
		# subject to
		# b = Ax
		# h >= Gx TODO: check this
		#
		# where b = 0

		## Create A matrix (stoichiometry + other things)

		self._A = cvxopt.spmatrix(values, rowIndexes, colIndexes)

		## Create objective function f

		objectiveFunction = np.zeros(self._nEdges, np.float64)

		objectiveFunction[objIndexes] = objValues

		self._f = cvxopt.matrix(-objectiveFunction) # negative, since GLPK minimizes

		self._lowerBound = np.zeros(self._nEdges, np.float64)
		self._upperBound = np.empty(self._nEdges, np.float64)

		self._upperBound.fill(np.inf)
		# TODO: something analogous for the lower bound

		self._G = cvxopt.matrix(np.concatenate(
			[np.identity(self._nEdges, np.float64), -np.identity(self._nEdges, np.float64)], axis = 0
			))

		self._b = cvxopt.matrix(np.zeros(self._nNodes, np.float64))

		## Create matrix for computing output
		self._outputCalcMatrix = -np.array(cvxopt.matrix(
			self._A[outputMoleculeIndexes, outputReactionIndexes]
			))


	def _edgeAdd(self, edgeName):
		if edgeName in self._edgeNames:
			raise FBAException("Edge already exists: {}".format(edgeName))

		else:
			self._edgeNames.append(edgeName)
			return len(self._edgeNames) - 1


	def _edgeIndex(self, edgeName, createIfDoesNotExists = False):
		try:
			return self._edgeNames.index(edgeName)

		except ValueError:
			if createIfDoesNotExists:
				return self._edgeAdd(edgeName)

			else:
				raise FBAException("Edge does not exist: {}".format(edgeName))


	def _nodeAdd(self, nodeName):
		if nodeName in self._nodeNames:
			raise FBAException("Node already exists: {}".format(nodeName))

		else:
			self._nodeNames.append(nodeName)
			return len(self._nodeNames) - 1


	def _nodeIndex(self, nodeName, createIfDoesNotExists = False):
		try:
			return self._nodeNames.index(nodeName)

		except ValueError:
			if createIfDoesNotExists:
				return self._nodeAdd(nodeName)

			else:
				raise FBAException("Node does not exist: {}".format(nodeName))


	# Constraint setup

	def externalMoleculeCountsIs(self, counts):
		self._lowerBound[self._externalExchangeIndexes] = -counts


	def internalMoleculeCountsIs(self, counts):
		raise NotImplementedError()


	def enzymeCountsIs(self, counts):
		raise NotImplementedError()


	# Evaluation

	def run(self):
		h = cvxopt.matrix(
			np.concatenate([self._upperBound, -self._lowerBound], axis = 0)
			)

		oldOptions = cvxopt.solvers.options.copy()

		cvxopt.solvers.options["LPX_K_MSGLEV"] = 0

		solution = cvxopt.solvers.lp(self._f, self._G, h, self._A, self._b, solver = "glpk")

		self._rawSolution = solution

		# TODO: raise/return flag on failed optimization

		self._edgeFluxes = np.array(self._rawSolution["x"]).flatten()


	# Output

	def outputMoleculeIDs(self):
		return self._outputMoleculeIDs


	def outputMoleculeCounts(self):
		# Must compute and return two (potentially overlapping) sets of 
		# molecule counts:
		# - internal input usage (TODO)
		# - objective output production

		# TODO: something about noninteger returns

		return np.dot(self._outputCalcMatrix, self._edgeFluxes[self._outputReactionIndexes])


if __name__ == "__main__":
	reactionData = {
		"A to B":{
			"enzymeID":None,
			"isReversible":False,
			"stoichiometry":{
				"A":-1,
				"B":+1,
				}
			},
		"AB2 to C":{
			"enzymeID":None,
			"isReversible":True,
			"stoichiometry":{
				"A":-1,
				"B":-2,
				"C":+1,
				}
			},
		}

	# TODO: split up stoichiometry/enzyme association/reversibility into separate arguments

	externalExchangedMolecules = ["A"]

	objective = {"B":20, "C":10}

	fba = FluxBalanceAnalysis(reactionData, externalExchangedMolecules, objective)

	fba.externalMoleculeCountsIs(10)
	fba.run()

	print fba.outputMoleculeIDs()
	print fba.outputMoleculeCounts()

