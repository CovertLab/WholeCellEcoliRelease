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

from collections import defaultdict

# import numpy as np
import cvxopt

class FBAException(Exception):
	pass


class AlreadyExistsException(FBAException):
	pass


class DoesNotExistException(FBAException):
	pass


class InvalidBoundaryException(FBAException):
	pass


class FluxBalanceAnalysis(object):
	""" FluxBalanceAnalysis

	Solver for various FBA implementations.

	
	Required arguments:

	- reactionStoich, a dict of strings:dicts (reactionID:reaction stoich)
		Each value in the dict is a dict of molecule ID to stoichiometry pairs

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

	- reversibleReactions, a list of strings (reactionIDs for reversible reactions)

	- reactionEnzymes, a dict of strings:strings (reactionID:enzymeID)
	
	- reactionRates, a dict of strings:floats (reactionID:catalytic rate constant * dt)
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

	# Format strings for autogenerated IDs

	## Reverse reactions
	_generatedID_reverseReaction = "{} (reverse)"

	## Exchange fluxes
	_generatedID_externalExchange = "{} external exchange"
	_generatedID_internalExchange = "{} internal exchange"

	## Objective
	_generatedID_moleculesToEquivalents = "molecules of {} to fractional objective equivalents"
	_generatedID_moleculeEquivalents = "fractional objective equivalent for {}"

	## Rate-constrained enzymes
	_generatedID_enzymeEquivRateConstrained = "{} equivalent (rate-constrained)"
	_generatedID_enzymeUsageRateConstrained = "{} usage (rate-constrained)"

	## Bool-constrained enzymes
	_generatedID_enzymeEquivBoolConstrained = "{} equivalent (bool-constrained)"
	_generatedID_enzymeUsageBoolConstrained = "{} usage (bool-constrained)"

	## Flex FBA
	_generatedID_fractionalDifferenceLeading = "difference between fractional objective equivalents of leading molecule and {}"
	_generatedID_fractionalDifferenceBiomass = "difference between fractional objective equivalents of {} and biomass objective"

	_generatedID_fractionsOut = "fractional objective equivalents of {} out"
	_generatedID_fractionalDifferenceLeadingOut = "difference between fractional objective equivalents of leading molecule and {} out"
	_generatedID_fractionalDifferenceBiomassOut = "difference between fractional objective equivalents of {} and biomass objective out"

	# Default values, for clarity
	_lowerBoundDefault = 0
	_upperBoundDefault = np.inf

	_standardObjectiveReactionName = "Standard biomass objective reaction"


	# Initialization

	def __init__(self, reactionStoich, externalExchangedMolecules, objective,
			objectiveType = None, objectiveParameters = None,
			internalExchangedMolecules = None, reversibleReactions = None, 
			reactionEnzymes = None, reactionRates = None, moleculeMasses = None):

		# Set up attributes

		self._nodeNames = []
		self._edgeNames = []

		# Set up running values for initialization

		## Used for creating the sparse matrix
		self._nodeIndexs = []
		self._edgeIndexes = []
		self._values = []

		## Used for objective and bounds
		self._objIndexes = []
		self._objValues = []

		self._lowerBoundIndexes = []
		self._lowerBoundValues = []

		self._upperBoundIndexes = []
		self._upperBoundValues = []

		## Output calculations
		self._outputMoleculeIndexes = []
		self._outputReactionIndexes = []

		# Set up reversible reactions
		if reversibleReactions is not None:
			for reactionID in reversibleReactions:
				reverseReactionID = self._generatedID_reverseReaction.format(reactionID)
				
				reactionStoich[reverseReactionID] = {
					moleculeID:-stoichCoeff
					for moleculeID, stoichCoeff in reactionStoich[reactionID].viewitems()
					}

				if reactionEnzymes is not None and reactionEnzymes.has_key(reactionID):
					reactionEnzymes[reverseReactionID] = reactionEnzymes[reactionID]

				if reactionRates is not None and reactionRates.has_key(reactionID):
					reactionRates[reverseReactionID] = reactionRates[reactionID]

		# Call indivdual initialization methods

		self._initReactionNetwork(reactionStoich)
		self._initExternalExchange(externalExchangedMolecules)

		self._initObjectiveEquivalents(objective)

		if objectiveType is None or objectiveType == "standard":
			self._initObjectiveStandard(objective)

		elif objectiveType == "flexible":
			self._initObjectiveFlexible(objective, objectiveParameters)

		elif objectiveType == "pools":
			self._initObjectivePools(objective)

		else:
			raise FBAException("Unrecognized objectiveType: {}".format(objectiveType))

		self._initInternalExchange(internalExchangedMolecules)

		self._initEnzymeConstraints(reactionEnzymes, reactionRates)

		self._initMass(moleculeMasses)

		# Finalize

		self._finalizeMisc()
		self._finalizeMatrices()

		# Set up values that will change between runs

		self.externalMoleculeLevelsIs(0)
		self.internalMoleculeLevelsIs(0)
		self.enzymeLevelsIs(0)


	def _initReactionNetwork(self, reactionStoich):
		""" Create the reaction network, initializing molecules and biochemical
		reactions. """
		
		reactionIndexes = []

		for reactionID, stoichiometry in reactionStoich.viewitems():
			colIndex = self._edgeAdd(reactionID)

			for moleculeID, stoichCoeff in stoichiometry.viewitems():
				rowIndex = self._nodeIndex(moleculeID, True)

				self._nodeIndexs.append(rowIndex)
				self._edgeIndexes.append(colIndex)
				self._values.append(stoichCoeff)

			reactionIndexes.append(colIndex)

		self._reactionIndexes = np.array(reactionIndexes)


	def _initExternalExchange(self, externalExchangedMolecules):
		"""Create external (media) exchange reactions."""

		externalExchangeIndexes = []
		externalMoleculeIDs = []

		for moleculeID in externalExchangedMolecules:
			exchangeID = self._generatedID_externalExchange.format(moleculeID)

			colIndex = self._edgeAdd(exchangeID)
			rowIndex = self._nodeIndex(moleculeID)

			# NOTE: The convention, if I am not mistaken, is to define 
			# exchanges as pointing out of the system, even though it is common
			# to think of exchanges as inputs.  Regardless this choice should 
			# have no impact outside of this class.

			self._nodeIndexs.append(rowIndex)
			self._edgeIndexes.append(colIndex)
			self._values.append(-1)

			externalExchangeIndexes.append(colIndex)
			externalMoleculeIDs.append(moleculeID)

		self._externalExchangeIndexes = np.array(externalExchangeIndexes, np.int64)
		self._externalMoleculeIDs = tuple(externalMoleculeIDs)


	def _initObjectiveEquivalents(self, objective):
		"""Create pseudo-reactions that convert molecules into their fractional
		objective equivalents.  The objectiveType determines how these 
		fractions are used."""

		objectiveEquivalentIndexes = []

		for moleculeID, coeff in objective.viewitems():
			molecule_rowIndex = self._nodeIndex(moleculeID)

			pseudoFluxID = self._generatedID_moleculesToEquivalents.format(moleculeID)
			colIndex = self._edgeAdd(pseudoFluxID)

			objectiveEquivID = self._generatedID_moleculeEquivalents.format(moleculeID)
			objectiveEquiv_rowIndex = self._nodeAdd(objectiveEquivID)

			self._nodeIndexs.append(molecule_rowIndex)
			self._edgeIndexes.append(colIndex)
			self._values.append(-coeff)

			self._nodeIndexs.append(objectiveEquiv_rowIndex)
			self._edgeIndexes.append(colIndex)
			self._values.append(+1)

			self._outputMoleculeIndexes.append(molecule_rowIndex)
			objectiveEquivalentIndexes.append(objectiveEquiv_rowIndex)

			self._outputReactionIndexes.append(colIndex)


	def _initObjectiveStandard(self, objective):
		"""Create the pseudo-reaction for the standard biomass objective.  In 
		the standard objective, all molecules must be created/destroyed in 
		prescribed ratios."""

		colIndex = self._edgeAdd(self._standardObjectiveReactionName)

		nObjectiveEquivalents = len(objectiveEquivalentIndexes)

		for moleculeID in objective.viewkeys():
			objectiveEquivID = self._generatedID_moleculeEquivalents.format(moleculeID)
			objectiveEquiv_rowIndex = self._nodeAdd(objectiveEquivID)

			self._nodeIndexs.append(objectiveEquiv_rowIndex)
			self._edgeIndexes.append(colIndex)
			self._values.append(-1)

		self._objIndexes.append(colIndex)
		self._objValues.append(+1)


	def _initObjectiveFlexible(self, objective, objectiveParameters):
		"""Create the abstractions needed for the flexFBA objective.  In brief,
		flexFBA permits partial biomass objective satisfaction for individual
		molecules if network disruptions inhibit molecule production."""

		# Load parameters
		leadingMoleculeID = objectiveParameters["leading molecule ID"]

		if not objective.has_key(leadingMoleculeID):
			raise FBAException("flexFBA leading molecule must be in the objective")

		fractionalDifferenceWeight = objectiveParameters["gamma"]

		if fractionalDifferenceWeight < 0:
			raise FBAException("flexFBA gamma paramter must be nonnegative")

		biomassSatisfactionWeight = objectiveParameters["beta"]

		if biomassSatisfactionWeight < 0:
			raise FBAException("flexFBA beta paramter must be nonnegative")

		biomass_colIndex = self._edgeAdd(self._standardObjectiveReactionName)

		# Add biomass to objective
		self._objIndexes.append(biomass_colIndex)
		self._objValues.append(biomassSatisfactionWeight)

		# Create fraction and biomass outputs
		for moleculeID in objective.viewkeys():
			fractionID = self._generatedID_moleculeEquivalents.format(moleculeID)
			fraction_rowIndex = self._nodeIndex(fractionID)

			# Biomass out
			self._nodeIndexs.append(fraction_rowIndex)
			self._edgeIndexes.append(biomass_colIndex)
			self._values.append(-1)

			# Fraction out
			fractionOutID = self._generatedID_fractionsOut.format(moleculeID)
			fractionOut_colIndex = self._edgeAdd(fractionOutID)

			self._nodeIndexs.append(fraction_rowIndex)
			self._edgeIndexes.append(fractionOut_colIndex)
			self._values.append(-1)

			if moleculeID == leadingMoleculeID:
				# Add leading molecule to objective
				self._objIndexes.append(fractionOut_colIndex)
				self._objValues.append(+1)

		# Create fraction differences (leading - other), used in objective and constraints
		leadingMoleculeToFractionID = self._generatedID_moleculesToEquivalents.format(leadingMoleculeID)
		leadingMoleculeToFraction_colIndex = self._edgeIndex(leadingMoleculeToFractionID)

		for moleculeID in objective.viewkeys():
			if moleculeID == leadingMoleculeID:
				continue

			fractionDifferenceLeadingID = self._generatedID_fractionalDifferenceLeading.format(moleculeID)
			fractionDifferenceLeading_rowIndex = self._nodeAdd(fractionDifferenceLeadingID)

			self._nodeIndexs.append(fractionDifferenceLeading_rowIndex)
			self._edgeIndexes.append(leadingMoleculeToFraction_colIndex)
			self._values.append(+1)
			
			moleculeToFractionID = self._generatedID_moleculesToEquivalents.format(moleculeID)
			moleculeToFraction_colIndex = self._edgeIndex(moleculeToFractionID)

			self._nodeIndexs.append(fractionDifferenceLeading_rowIndex)
			self._edgeIndexes.append(moleculeToFraction_colIndex)
			self._values.append(-1)

			fractionDifferenceLeadingOutID = self._generatedID_fractionalDifferenceLeadingOut.format(moleculeID)
			fractionDifferenceLeadingOut_colIndex = self._edgeAdd(fractionDifferenceLeadingOutID)

			self._nodeIndexs.append(fractionDifferenceLeading_rowIndex)
			self._edgeIndexes.append(fractionDifferenceLeadingOut_colIndex)
			self._values.append(-1)

			self._objIndexes.append(fractionDifferenceLeadingOut_colIndex)
			self._objValues.append(-fractionalDifferenceWeight)

		# Create biomass differences (fraction - biomass), used in constraints

		for moleculeID in objective.viewkeys():
			fractionDifferenceBiomassID = self._generatedID_fractionalDifferenceBiomass.format(moleculeID)
			fractionDifferenceBiomass_rowIndex = self._nodeAdd(fractionDifferenceBiomassID)

			moleculeToFractionID = self._generatedID_moleculesToEquivalents.format(moleculeID)
			moleculeToFraction_colIndex = self._edgeIndex(moleculeToFractionID)

			self._nodeIndexs.append(fractionDifferenceBiomass_rowIndex)
			self._edgeIndexes.append(moleculeToFraction_colIndex)
			self._values.append(+1)

			self._nodeIndexs.append(fractionDifferenceBiomass_rowIndex)
			self._edgeIndexes.append(biomass_colIndex)
			self._values.append(-1)

			fractionDifferenceBiomassOutID = self._generatedID_fractionalDifferenceBiomassOut.format(moleculeID)
			fractionDifferenceBiomassOut_colIndex = self._edgeAdd(fractionDifferenceBiomassOutID)

			self._nodeIndexs.append(fractionDifferenceBiomass_rowIndex)
			self._edgeIndexes.append(fractionDifferenceBiomassOut_colIndex)
			self._values.append(-1)


	def _initObjectivePools(self, objective):
		"""Not implemented"""
		raise NotImplementedError()


	def _initInternalExchange(self, internalExchangedMolecules):
		"""Create internal (byproduct) exchange reactions."""

		internalExchangeIndexes = []
		internalMoleculeIDs = []

		if internalExchangedMolecules is not None:
			for moleculeID in internalExchangedMolecules:
				exchangeID = self._generatedID_internalExchange.format(moleculeID)

				colIndex = self._edgeAdd(exchangeID)
				rowIndex = self._nodeIndex(moleculeID)

				self._nodeIndexs.append(rowIndex)
				self._edgeIndexes.append(colIndex)
				self._values.append(-1)

				internalExchangeIndexes.append(colIndex)
				internalMoleculeIDs.append(moleculeID)

				self._outputMoleculeIndexes.append(rowIndex)
				self._outputReactionIndexes.append(colIndex)

		self._internalExchangeIndexes = np.array(internalExchangeIndexes, np.int64)
		self._internalMoleculeIDs = tuple(internalMoleculeIDs)


	def _initEnzymeConstraints(self, reactionEnzymes, reactionRates):
		"""Create abstractions needed to constrain metabolic reactions by 
		enzyme availability.

		There are two types of enzyme restrictions.  If the catalytic rate is 
		unknown, any non-zero level of the enzyme allows the reaction to
		proceed.  This is a "boolean" constraint.  If the catalytic rate is
		known, then reactions may proceed up to their kinetic limit.  Reactions
		can share enzymes.  There is currently no support for reactions with 
		multiple annotated enzymes."""

		enzymeUsageRateConstrainedIndexes = []
		enzymeUsageBoolConstrainedIndexes = []

		if reactionEnzymes is not None:

			## First create the pseudometabolites and enzyme usage columns
			self._enzymeIDs = tuple(set(reactionEnzymes.values()))

			for enzymeID in self._enzymeIDs:
				# Create pseudometabolite and flux for rate-constrained
				enzymeEquivalentRateID = self._generatedID_enzymeEquivRateConstrained.format(enzymeID)
				enzymeEquivalentRate_rowIndex = self._nodeAdd(enzymeEquivalentRateID)

				enzymeUsageRateID = self._generatedID_enzymeUsageRateConstrained.format(enzymeID)
				enzymeUsageRate_colIndex = self._edgeAdd(enzymeUsageRateID)

				self._nodeIndexs.append(enzymeEquivalentRate_rowIndex)
				self._edgeIndexes.append(enzymeUsageRate_colIndex)
				self._values.append(+1)

				enzymeUsageRateConstrainedIndexes.append(enzymeUsageRate_colIndex)

				# Create pseudometabolite and flux for bool-constrained
				enzymeEquivalentBoolID = self._generatedID_enzymeEquivBoolConstrained.format(enzymeID)
				enzymeEquivalentBool_rowIndex = self._nodeAdd(enzymeEquivalentBoolID)

				enzymeUsageBoolID = self._generatedID_enzymeUsageBoolConstrained.format(enzymeID)
				enzymeUsageBool_colIndex = self._edgeAdd(enzymeUsageBoolID)

				self._nodeIndexs.append(enzymeEquivalentBool_rowIndex)
				self._edgeIndexes.append(enzymeUsageBool_colIndex)
				self._values.append(+1)

				enzymeUsageBoolConstrainedIndexes.append(enzymeUsageBool_colIndex)


			for reactionID, enzymeID in reactionEnzymes.viewitems():
				reaction_colIndex = self._edgeIndex(reactionID)

				if reactionRates is not None and reactionRates.has_key(reactionID):
					reactionRate = reactionRates[reactionID]

					if reactionRate <= 0:
						raise FBAException("Reaction rates must be positive ({}, {})".format(reactionID, reactionRate))

					enzymesPerReaction = -1/reactionRate
					enzymeEquivalentRateID = self._generatedID_enzymeEquivRateConstrained.format(enzymeID)
					enzymeEquivalentRate_rowIndex = self._nodeIndex(enzymeEquivalentRateID)

					self._nodeIndexs.append(enzymeEquivalentRate_rowIndex)
					self._edgeIndexes.append(reaction_colIndex)
					self._values.append(enzymesPerReaction)
					

				else:
					enzymeEquivalentBoolID = self._generatedID_enzymeEquivBoolConstrained.format(enzymeID)
					enzymeEquivalentBool_rowIndex = self._nodeIndex(enzymeEquivalentBoolID)

					self._nodeIndexs.append(enzymeEquivalentBool_rowIndex)
					self._edgeIndexes.append(reaction_colIndex)
					self._values.append(-1)

		self._enzymeUsageRateConstrainedIndexes = np.array(enzymeUsageRateConstrainedIndexes, np.int64)
		self._enzymeUsageBoolConstrainedIndexes = np.array(enzymeUsageBoolConstrainedIndexes, np.int64)


	def _initMass(self, moleculeMasses):
		"""Not implemented.

		Tracking the mass entering the system through metabolism is crucial for
		insuring closure. It can also be used to constrain the maximal rate of
		growth."""

		if moleculeMasses is not None:
			raise NotImplementedError()


	def _finalizeMisc(self):
		"""Finalize values accumulated during the various _init* methods."""
		self._nEdges = len(self._edgeNames)
		self._nNodes = len(self._nodeNames)

		self._outputMoleculeIDs = tuple([self._nodeNames[index] for index in self._outputMoleculeIndexes])

		self._outputReactionIndexes = np.array(self._outputReactionIndexes)


	def _finalizeMatrices(self):
		"""Create the abstractions needed for linear programming and for 
		computing various outputs.

		Optimization problem:
		
		\max_v f^T x
		
		subject to
		b = Ax
		h >= Gx
		
		where b = 0"""

		# Create A matrix (stoichiometry + other things)

		self._A = cvxopt.spmatrix(self._values, self._nodeIndexs, self._edgeIndexes)

		# Create objective function f

		objectiveFunction = np.zeros(self._nEdges, np.float64)

		objectiveFunction[self._objIndexes] = self._objValues

		self._f = cvxopt.matrix(-objectiveFunction) # negative, since GLPK minimizes (TODO: max/min mode setting)

		self._lowerBound = np.empty(self._nEdges, np.float64)
		self._lowerBound.fill(self._lowerBoundDefault)
		self._lowerBound[self._lowerBoundIndexes] = self._lowerBoundValues

		self._upperBound = np.empty(self._nEdges, np.float64)
		self._upperBound.fill(self._upperBoundDefault)
		self._upperBound[self._upperBoundIndexes] = self._upperBoundValues

		self._G = cvxopt.matrix(np.concatenate(
			[np.identity(self._nEdges, np.float64), -np.identity(self._nEdges, np.float64)], axis = 0
			))

		self._b = cvxopt.matrix(np.zeros(self._nNodes, np.float64))

		# Create matrix for computing output
		self._outputCalcMatrix = -np.array(cvxopt.matrix(self._A[
			self._outputMoleculeIndexes, self._outputReactionIndexes.tolist() # NOTE: this type of indexing REQUIRES lists, not arrays
			]))


	def _edgeAdd(self, edgeName):
		if edgeName in self._edgeNames:
			raise AlreadyExistsException("Edge already exists: {}".format(edgeName))

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
				raise DoesNotExistException("Edge does not exist: {}".format(edgeName))


	def _nodeAdd(self, nodeName):
		if nodeName in self._nodeNames:
			raise AlreadyExistsException("Node already exists: {}".format(nodeName))

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
				raise DoesNotExistException("Node does not exist: {}".format(nodeName))


	# Constraint setup

	# NOTE: I've used the nonspecific term "levels" since these could be
	# mole-concentration, number-concentration, or counts, depending on the 
	# formulation of the problem.  All that matters is that initialization 
	# parameters have consistent units.

	def externalMoleculeIDs(self):
		return self._externalMoleculeIDs


	def externalMoleculeLevelsIs(self, levels):
		levels = np.array(levels)
		if (levels < 0).any():
			raise InvalidBoundaryException("Negative molecule levels not allowed")

		self._lowerBound[self._externalExchangeIndexes] = -levels


	def internalMoleculeIDs(self):
		return self._internalMoleculeIDs


	def internalMoleculeLevelsIs(self, levels):
		levels = np.array(levels)
		if (levels < 0).any():
			raise InvalidBoundaryException("Negative molecule levels not allowed")

		self._lowerBound[self._internalExchangeIndexes] = -levels


	def enzymeIDs(self):
		return self._enzymeIDs


	def enzymeLevelsIs(self, levels):
		if self._enzymeUsageRateConstrainedIndexes.size == 0:
			return

		levels = np.array(levels)
		if (levels < 0).any():
			raise InvalidBoundaryException("Negative enzyme levels not allowed")

		# Rate-constrained
		self._upperBound[self._enzymeUsageRateConstrainedIndexes] = levels

		# Boolean-constrained (enzyme w/o an annotated rate)
		boolConstraint = np.zeros(self._enzymeUsageBoolConstrainedIndexes.size, np.float64)
		boolConstraint[levels > 0] = np.inf
		self._upperBound[self._enzymeUsageBoolConstrainedIndexes] = boolConstraint


	def maxReactionFluxIs(self, reactionID, maxFlux, raiseForReversible = True):
		colIndex = self._edgeIndex(reactionID)

		if maxFlux < 0:
			raise InvalidBoundaryException("Maximum reaction flux must be at least 0")

		if maxFlux < self._lowerBound[colIndex]:
			raise InvalidBoundaryException("Maximum reaction flux must be greater than or equal to the minimum flux")

		reverseReactionID = self._generatedID_reverseReaction.format(reactionID)

		if raiseForReversible and reverseReactionID in self._edgeNames:
			raise FBAException((
				"Setting the maximum reaction flux is ambiguous since " +
				"reaction {} has both a forward [{}] and reverse [{}] " +
				"component.  Call this method with argument " + 
				"raiseForReversible = False if this is intended behavior."
				).format(reactionID, reactionID, reverseReactionID))

		self._upperBound[colIndex] = maxFlux


	def minReactionFluxIs(self, reactionID, minFlux, raiseForReversible = True):
		colIndex = self._edgeIndex(reactionID)

		if minFlux < 0:
			raise InvalidBoundaryException("Maximum reaction flux must be at least 0")

		if minFlux > self._upperBound[colIndex]:
			raise InvalidBoundaryException("Minimum reaction flux must be less than or equal to the maximum flux")

		reverseReactionID = self._generatedID_reverseReaction.format(reactionID)

		if raiseForReversible and reverseReactionID in self._edgeNames:
			raise FBAException((
				"Setting the minimum reaction flux is ambiguous since " +
				"reaction {} has both a forward [{}] and reverse [{}] " +
				"component.  Call this method with argument " + 
				"raiseForReversible = False if this is intended behavior."
				).format(reactionID, reactionID, reverseReactionID))

		self._lowerBound[colIndex] = minFlux


	# Evaluation

	def run(self):
		h = cvxopt.matrix(
			np.concatenate([self._upperBound, -self._lowerBound], axis = 0)
			)

		oldOptions = cvxopt.solvers.options.copy()

		cvxopt.solvers.options["LPX_K_MSGLEV"] = 0

		solution = cvxopt.solvers.lp(self._f, self._G, h, self._A, self._b, solver = "glpk")

		cvxopt.solvers.options.update(oldOptions)

		self._rawSolution = solution

		# TODO: raise/return flag on failed optimization

		self._edgeFluxes = np.array(self._rawSolution["x"]).flatten()


	# Output

	def outputMoleculeIDs(self):
		return self._outputMoleculeIDs


	def outputMoleculeLevelsChange(self):
		# Must compute and return two (potentially overlapping) sets of 
		# molecule counts:
		# - internal input usage (TODO)
		# - objective output production

		return np.dot(self._outputCalcMatrix, self._edgeFluxes[self._outputReactionIndexes])


	def externalExchangeFlux(self, moleculeID):
		return -self._edgeFluxes[
			self._edgeIndex(self._generatedID_externalExchange.format(moleculeID))
			]


	def internalExchangeFlux(self, moleculeID):
		return -self._edgeFluxes[
			self._edgeIndex(self._generatedID_internalExchange.format(moleculeID))
			]


	def reactionFlux(self, reactionID):
		return self._edgeFluxes[self._edgeIndex(reactionID)]


	def objectiveReactionFlux(self):
		try:
			colIndex = self._edgeIndex(self._standardObjectiveReactionName)

		except DoesNotExistException:
			raise FBAException("No objective reaction flux implemented for this solver type")

		return self._edgeFluxes[colIndex]


# Test data

def loadKB():
	from wholecell.reconstruction.knowledge_base_ecoli import KnowledgeBaseEcoli
	return KnowledgeBaseEcoli()


def setupFeist(kb):
	# Create the KB and parse into FBA inputs

	objectiveRaw = {
		'10fthf[c]' : -0.000223,
		'2ohph[c]' : -0.000223,
		'adp[c]' : 59.810000000000002,
		'ala-L[c]' : -0.51370000000000005,
		'amet[c]' : -0.000223,
		'arg-L[c]' : -0.29580000000000001,
		'asn-L[c]' : -0.24110000000000001,
		'asp-L[c]' : -0.24110000000000001,
		'atp[c]' : -59.984000000000002,
		'ca2[c]' : -0.0047369999999999999,
		'cl[c]' : -0.0047369999999999999,
		'coa[c]' : -0.00057600000000000001,
		'cobalt2[c]' : -0.0031580000000000002,
		'ctp[c]' : -0.13350000000000001,
		'cu2[c]' : -0.0031580000000000002,
		'cys-L[c]' : -0.091579999999999995,
		'datp[c]' : -0.026169999999999999,
		'dctp[c]' : -0.027019999999999999,
		'dgtp[c]' : -0.027019999999999999,
		'dttp[c]' : -0.026169999999999999,
		'fad[c]' : -0.000223,
		'fe2[c]' : -0.0071060000000000003,
		'fe3[c]' : -0.0071060000000000003,
		'gln-L[c]' : -0.26319999999999999,
		'glu-L[c]' : -0.26319999999999999,
		'gly[c]' : -0.61260000000000003,
		'gtp[c]' : -0.21510000000000001,
		'h2o[c]' : -54.462000000000003,
		'h[c]' : 59.810000000000002,
		'his-L[c]' : -0.094740000000000005,
		'ile-L[c]' : -0.29049999999999998,
		'k[c]' : -0.17760000000000001,
		'kdo2lipid4[e]' : -0.019449999999999999,
		'leu-L[c]' : -0.45050000000000001,
		'lys-L[c]' : -0.34320000000000001,
		'met-L[c]' : -0.1537,
		'mg2[c]' : -0.0078949999999999992,
		'mlthf[c]' : -0.000223,
		'mn2[c]' : -0.0031580000000000002,
		'mobd[c]' : -0.0031580000000000002,
		'murein5px4p[p]' : -0.01389,
		'nad[c]' : -0.0018309999999999999,
		'nadp[c]' : -0.00044700000000000002,
		'nh4[c]' : -0.011842999999999999,
		'pe160[c]' : -0.022329999999999999,
		'pe160[p]' : -0.041480000000000003,
		'pe161[c]' : -0.02632,
		'pe161[p]' : -0.048890000000000003,
		'phe-L[c]' : -0.1759,
		'pheme[c]' : -0.000223,
		'pi[c]' : 59.805999999999997,
		'ppi[c]' : 0.77390000000000003,
		'pro-L[c]' : -0.22109999999999999,
		'pydx5p[c]' : -0.000223,
		'ribflv[c]' : -0.000223,
		'ser-L[c]' : -0.21579999999999999,
		'sheme[c]' : -0.000223,
		'so4[c]' : -0.0039480000000000001,
		'thf[c]' : -0.000223,
		'thmpp[c]' : -0.000223,
		'thr-L[c]' : -0.25369999999999998,
		'trp-L[c]' : -0.056840000000000002,
		'tyr-L[c]' : -0.13789999999999999,
		'udcpdp[c]' : -5.5000000000000002e-05,
		'utp[c]' : -0.14410000000000001,
		'val-L[c]' : -0.42320000000000002,
		'zn2[c]' : -0.0031580000000000002,
		}

	objective = {}
	for moleculeID_raw, coeff in objectiveRaw.viewitems():
		moleculeID = moleculeID_raw[:-2].upper() + moleculeID_raw[-2:]

		if moleculeID == "KDO2LIPID4[e]":
			moleculeID = "KDO2LIPID4[o]"

		objective[moleculeID] = -coeff

	import re
	rxns = [
		x for x in kb.metabolismBiochemicalReactions
		if (
			not re.match(".*_[0-9]$", x["id"])
			or x["id"].endswith("_0")
			or "PFK_2" in x["id"]
			)
		]

	reactionStoich = {
		rxn["id"]:
		{entry["molecule"]+"["+entry["location"]+"]":entry["coeff"] for entry in rxn["stoichiometry"]}
		for rxn in rxns
		if len(rxn["stoichiometry"]) > 1 # no exchange reactions!
		}

	reversibleReactions = [rxn["id"] for rxn in rxns if not rxn["dir"]]

	mediaEx = kb.metabolismMediaEx

	externalExchangedMolecules = [rxn["met"] for rxn in mediaEx]

	atpId = "ATP[c]"

	# Create FBA instance
	fba = FluxBalanceAnalysis(
		reactionStoich,
		externalExchangedMolecules,
		objective,
		objectiveType = "flexible",
		objectiveParameters = {
			"gamma":0.1,
			"beta":100,
			"leading molecule ID":atpId
			},
		reversibleReactions = reversibleReactions,
		)

	# Set constraints
	## External molecules
	externalMoleculeIDs = fba.externalMoleculeIDs()

	unconstrainedExchange = ("CA2[e]", "CL[e]", "CO2[e]", "COBALT2[e]", 
		"CU2[e]", "FE2[e]", "FE3[e]", "H[e]", "H2O[e]", "K[e]", "MG2[e]",
		"MN2[e]", "MOBD[e]", "NA1[e]", "NH4[e]", "PI[e]", "SO4[e]", "TUNGS[e]",
		"ZN2[e]")

	constrainedExchange = {
		"CBL1[e]":0.01,
		"GLC-D[e]":8,
		"O2[e]":18.5
		}

	externalMoleculeLevels = np.zeros(len(externalMoleculeIDs), np.float64)

	for index, moleculeID in enumerate(externalMoleculeIDs):
		if moleculeID in unconstrainedExchange:
			externalMoleculeLevels[index] = np.inf

		elif constrainedExchange.has_key(moleculeID):
			externalMoleculeLevels[index] = constrainedExchange[moleculeID]

	fba.externalMoleculeLevelsIs(externalMoleculeLevels)

	## Set Feist's forced reactions

	### ATP maintenance
	fba.minReactionFluxIs("FEIST_ATPM", 8.39)
	fba.maxReactionFluxIs("FEIST_ATPM", 8.39)

	### Arbitrarily disabled reactions
	disabledReactions = ("FEIST_CAT_0", "FEIST_SPODM_0", "FEIST_SPODMpp", "FEIST_FHL_0_0", "FEIST_FHL_1_0")

	for reactionID in disabledReactions:
		fba.maxReactionFluxIs(reactionID, 0)

	return fba


def compareFeistToExpected():
	kb = loadKB()

	fba = setupFeist(kb)

	# Run model
	fba.run()

	# Output
	print fba.externalExchangeFlux("GLC-D[e]"), "(expected 8.)"
	print fba.externalExchangeFlux("O2[e]"), "(expected 16.27631182)"
	print fba.objectiveReactionFlux(), "(expected 0.73645239)"


if __name__ == "__main__":
	compareFeistToExpected()
