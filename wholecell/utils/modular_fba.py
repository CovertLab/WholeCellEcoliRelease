#!/usr/bin/env python

"""
@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 7/14/2014
"""

from __future__ import division

from collections import defaultdict
from itertools import izip
import warnings

import numpy as np

SOLVERS = {}
S_GUROBI = "gurobi"
S_GLPK = "glpk"
_SOLVER_PREFERENCE = (
	S_GUROBI,
	S_GLPK
	)

try:
	from ._netflow.nf_gurobi import NetworkFlowGurobi

except ImportError:
	pass

except Exception as e:
	# If this is a GurobiError, proceed without using gurobi, warning the user.
	if str(type(e)) == "<class 'gurobipy.GurobiError'>":
		print "GurobiError - gurobi will not be used."
	# Otherwise, raise the exception as normal
	else:
		raise e

else:
	SOLVERS[S_GUROBI] = NetworkFlowGurobi

try:
	from ._netflow.nf_glpk import NetworkFlowGLPK

except ImportError:
	pass

else:
	SOLVERS[S_GLPK] = NetworkFlowGLPK

if not SOLVERS:
	raise Exception("No solvers available.")

for solver in _SOLVER_PREFERENCE:
	if solver in SOLVERS:
		DEFAULT_SOLVER = solver
		break

else:
	raise Exception("Could not choose a default solver.")

# Errors

class FBAError(Exception):
	pass


class AlreadyExistsError(FBAError):
	pass


class DoesNotExistError(FBAError):
	pass


class InvalidBoundaryError(FBAError):
	pass


class SolverUnavailableError(FBAError):
	pass

# Classes


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

	## Pools FBA
	_generatedID_fractionBelowUnityOut = "fraction {} below unity, out"
	_generatedID_fractionAboveUnityOut = "fraction {} above unity, out"

	# Default values, for clarity
	_lowerBoundDefault = 0
	_upperBoundDefault = np.inf

	_standardObjectiveReactionName = "Standard biomass objective reaction"
	_massID = "Mass"
	_massOutName = "Mass out"

	_forcedUnityColName = "Column forced at unity"

	# Initialization

	def __init__(self, reactionStoich, externalExchangedMolecules, objective,
			objectiveType = None, objectiveParameters = None,
			internalExchangedMolecules = None, reversibleReactions = None,
			reactionEnzymes = None, reactionRates = None,
			moleculeMasses = None, maintenanceCost = None,
			maintenanceReaction = None,
			solver = DEFAULT_SOLVER):

		if solver not in SOLVERS:
			raise SolverUnavailableError(
				"Unrecognized or unavailable solver: {}".format(solver)
				)

		self._solver = SOLVERS[solver]()

		self._forceInternalExchange = False

		# Output calculations
		self._outputMoleculeIDs = []
		self._outputMoleculeCoeffs = []

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

			if internalExchangedMolecules is not None:
				raise FBAError(
					"Internal exchange molecules are automatically defined when using objectiveType = \"pools\""
					)

			internalExchangedMolecules = objective.keys()

		else:
			raise FBAError("Unrecognized objectiveType: {}".format(objectiveType))

		self._initInternalExchange(internalExchangedMolecules)

		self._initEnzymeConstraints(reactionEnzymes, reactionRates)

		self._initMass(externalExchangedMolecules, moleculeMasses)

		self._initMaintenance(maintenanceCost, maintenanceReaction)

		# Set up values that will change between runs

		self.externalMoleculeLevelsIs(0)
		self.internalMoleculeLevelsIs(0)
		self.enzymeLevelsIs(0)


	def _initReactionNetwork(self, reactionStoich):
		""" Create the reaction network, initializing molecules and biochemical
		reactions. """

		reactionIDs = []

		for reactionID, stoichiometry in reactionStoich.viewitems():
			for moleculeID, stoichCoeff in stoichiometry.viewitems():
				self._solver.flowMaterialCoeffIs(
					reactionID,
					moleculeID,
					stoichCoeff
					)

			reactionIDs.append(reactionID)

		self._reactionIDs = tuple(reactionIDs)


	def _initExternalExchange(self, externalExchangedMolecules):
		"""Create external (media) exchange reactions."""

		externalMoleculeIDs = []
		externalExchangeIDs = []

		for moleculeID in externalExchangedMolecules:
			exchangeID = self._generatedID_externalExchange.format(moleculeID)

			# NOTE: The convention, if I am not mistaken, is to define
			# exchanges as pointing out of the system, even though it is common
			# to think of exchanges as inputs.  Regardless this choice should
			# have no impact outside of this class.

			self._solver.flowMaterialCoeffIs(
				exchangeID,
				moleculeID,
				-1
				)

			externalMoleculeIDs.append(moleculeID)
			externalExchangeIDs.append(exchangeID)

		self._externalMoleculeIDs = tuple(externalMoleculeIDs)
		self._externalExchangeIDs = tuple(externalExchangeIDs)


	def _initObjectiveEquivalents(self, objective):
		"""Create pseudo-reactions that convert molecules into their fractional
		objective equivalents.  The objectiveType determines how these
		fractions are used."""

		for moleculeID, coeff in objective.viewitems():
			if coeff == 0:
				raise FBAError("Invalid objective coefficient - must be non-zero")

			pseudoFluxID = self._generatedID_moleculesToEquivalents.format(moleculeID)

			objectiveEquivID = self._generatedID_moleculeEquivalents.format(moleculeID)

			self._solver.flowMaterialCoeffIs(
				pseudoFluxID,
				moleculeID,
				-coeff
				)

			self._solver.flowMaterialCoeffIs(
				pseudoFluxID,
				objectiveEquivID,
				+1
				)

			# TODO: functionalize
			try:
				i = self._outputMoleculeIDs.index(moleculeID)

			except ValueError:
				self._outputMoleculeIDs.append(moleculeID)
				self._outputMoleculeCoeffs.append(dict())
				i = len(self._outputMoleculeIDs) - 1

			self._outputMoleculeCoeffs[i][pseudoFluxID] = -coeff


	def _initObjectiveStandard(self, objective):
		"""Create the pseudo-reaction for the standard biomass objective.  In
		the standard objective, all molecules must be created/destroyed in
		prescribed ratios."""

		for moleculeID in objective.viewkeys():
			objectiveEquivID = self._generatedID_moleculeEquivalents.format(moleculeID)

			self._solver.flowMaterialCoeffIs(
				self._standardObjectiveReactionName,
				objectiveEquivID,
				-1
				)

		self._solver.flowObjectiveCoeffIs(
			self._standardObjectiveReactionName,
			+1
			)


	def _initObjectiveFlexible(self, objective, objectiveParameters):
		"""Create the abstractions needed for the flexFBA objective.  In brief,
		flexFBA permits partial biomass objective satisfaction for individual
		molecules if network disruptions inhibit molecule production."""

		# Load parameters
		leadingMoleculeID = objectiveParameters["leading molecule ID"]

		if not objective.has_key(leadingMoleculeID):
			raise FBAError("flexFBA leading molecule must be in the objective")

		fractionalDifferenceWeight = objectiveParameters["gamma"]

		if fractionalDifferenceWeight < 0:
			raise FBAError("flexFBA gamma paramter must be nonnegative")

		biomassSatisfactionWeight = objectiveParameters["beta"]

		if biomassSatisfactionWeight < 0:
			raise FBAError("flexFBA beta paramter must be nonnegative")

		if any(coeff < 0 for coeff in objective.viewvalues()):
			warnings.warn("flexFBA is not designed to use negative biomass coefficients")

		# Add biomass to objective
		self._solver.flowObjectiveCoeffIs(
			self._standardObjectiveReactionName,
			biomassSatisfactionWeight
			)

		# Create fraction and biomass outputs
		for moleculeID in objective.viewkeys():
			fractionID = self._generatedID_moleculeEquivalents.format(moleculeID)

			# Biomass out
			self._solver.flowMaterialCoeffIs(
				self._standardObjectiveReactionName,
				fractionID,
				-1
				)

			# Fraction out
			fractionOutID = self._generatedID_fractionsOut.format(moleculeID)

			self._solver.flowMaterialCoeffIs(
				fractionOutID,
				fractionID,
				-1
				)

			if moleculeID == leadingMoleculeID:
				# Add leading molecule to objective
				self._solver.flowObjectiveCoeffIs(
					fractionOutID,
					+1
					)

		# Create fraction differences (leading - other), used in objective and constraints
		leadingMoleculeToFractionID = self._generatedID_moleculesToEquivalents.format(leadingMoleculeID)

		for moleculeID in objective.viewkeys():
			if moleculeID == leadingMoleculeID:
				continue

			fractionDifferenceLeadingID = self._generatedID_fractionalDifferenceLeading.format(moleculeID)

			self._solver.flowMaterialCoeffIs(
				leadingMoleculeToFractionID,
				fractionDifferenceLeadingID,
				+1
				)

			moleculeToFractionID = self._generatedID_moleculesToEquivalents.format(moleculeID)

			self._solver.flowMaterialCoeffIs(
				moleculeToFractionID,
				fractionDifferenceLeadingID,
				-1
				)

			fractionDifferenceLeadingOutID = self._generatedID_fractionalDifferenceLeadingOut.format(moleculeID)

			self._solver(
				fractionDifferenceLeadingOutID,
				fractionDifferenceLeadingID,
				-1
				)

			self._solver.flowObjectiveCoeffIs(
				fractionDifferenceLeadingOutID,
				-fractionalDifferenceWeight
				)

		# Create biomass differences (fraction - biomass), used in constraints

		for moleculeID in objective.viewkeys():
			fractionDifferenceBiomassID = self._generatedID_fractionalDifferenceBiomass.format(moleculeID)

			moleculeToFractionID = self._generatedID_moleculesToEquivalents.format(moleculeID)

			self._solver.flowMaterialCoeffIs(
				moleculeToFractionID,
				fractionDifferenceBiomassID,
				+1
				)

			self._solver.flowMaterialCoeffIs(
				fractionDifferenceBiomassID,
				self._standardObjectiveReactionName,
				-1
				)

			fractionDifferenceBiomassOutID = self._generatedID_fractionalDifferenceBiomassOut.format(moleculeID)

			self._solver.flowMaterialCoeffIs(
				fractionDifferenceBiomassOutID,
				fractionDifferenceBiomassID,
				-1
				)


	def _initObjectivePools(self, objective):
		"""Create the abstractions needed for FBA with pools.  The objective is
		to minimize the distance between the current metabolite level and some
		target level, as defined in the objective."""

		if any(coeff < 0 for coeff in objective.viewvalues()):
			raise FBAError("FBA with pools is not designed to use negative biomass coefficients")

		self._solver.maximizeObjective(False)
		self._forceInternalExchange = True

		# By forcing a column to be at unity, we can keep the definition of
		# the problem as b=Av where b=0.

		self._solver.flowLowerBoundIs(
			self._forcedUnityColName,
			+1
			)

		self._solver.flowUpperBoundIs(
			self._forcedUnityColName,
			+1
			)

		# Minimizing an absolute value requires splitting the term into two,
		# one for the positive values and one for the negative.

		for moleculeID in objective.viewkeys():
			objectiveEquivID = self._generatedID_moleculeEquivalents.format(moleculeID)

			# Add the forced -1 term so that we can define x_i = f_i - 1

			self._solver.flowMaterialCoeffIs(
				self._forcedUnityColName,
				objectiveEquivID,
				-1
				)

			# Add the term for when the flux out is below the expected value

			belowUnityID = self._generatedID_fractionBelowUnityOut.format(moleculeID)

			self._solver.flowMaterialCoeffIs(
				belowUnityID,
				objectiveEquivID,
				+1
				)

			self._solver.flowObjectiveCoeffIs(
				belowUnityID,
				+1
				)

			# Add the term for when the flux out is above the expected value

			aboveUnityID = self._generatedID_fractionAboveUnityOut.format(moleculeID)

			self._solver.flowMaterialCoeffIs(
				aboveUnityID,
				objectiveEquivID,
				-1
				)

			self._solver.flowObjectiveCoeffIs(
				aboveUnityID,
				+1
				)


	def _initInternalExchange(self, internalExchangedMolecules):
		"""Create internal (byproduct) exchange reactions."""

		internalMoleculeIDs = []

		if internalExchangedMolecules is not None:
			for moleculeID in internalExchangedMolecules:
				exchangeID = self._generatedID_internalExchange.format(moleculeID)

				self._solver.flowMaterialCoeffIs(
					exchangeID,
					moleculeID,
					-1
					)

				internalMoleculeIDs.append(moleculeID)

				# TODO: functionalize
				try:
					i = self._outputMoleculeIDs.index(moleculeID)

				except ValueError:
					self._outputMoleculeIDs.append(moleculeID)
					self._outputMoleculeCoeffs.append(dict())
					i = len(self._outputMoleculeIDs) - 1

				self._outputMoleculeCoeffs[i][exchangeID] = -1

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

		self._rateConstrainedEnzymeIDs = []

		if reactionEnzymes is not None:

			## First create the pseudometabolites and enzyme usage columns
			self._enzymeIDs = tuple(set(reactionEnzymes.values()))

			for enzymeID in self._enzymeIDs:
				# Create pseudometabolite and flux for rate-constrained
				enzymeEquivalentRateID = self._generatedID_enzymeEquivRateConstrained.format(enzymeID)
				enzymeUsageRateID = self._generatedID_enzymeUsageRateConstrained.format(enzymeID)

				self._solver.flowMaterialCoeffIs(
					enzymeUsageRateID,
					enzymeEquivalentRateID,
					+1
					)

				# Create pseudometabolite and flux for bool-constrained
				enzymeEquivalentBoolID = self._generatedID_enzymeEquivBoolConstrained.format(enzymeID)
				enzymeUsageBoolID = self._generatedID_enzymeUsageBoolConstrained.format(enzymeID)

				self._solver.flowMaterialCoeffIs(
					enzymeUsageBoolID,
					enzymeEquivalentBoolID,
					+1
					)


			for reactionID, enzymeID in reactionEnzymes.viewitems():
				if reactionRates is not None and reactionRates.has_key(reactionID):
					reactionRate = reactionRates[reactionID]

					if reactionRate <= 0:
						raise FBAError("Reaction rates must be positive ({}, {})".format(reactionID, reactionRate))

					enzymesPerReaction = -1/reactionRate
					enzymeEquivalentRateID = self._generatedID_enzymeEquivRateConstrained.format(enzymeID)

					self._solver.flowMaterialCoeffIs(
						reactionID,
						enzymeEquivalentRateID,
						enzymesPerReaction
						)

					self._rateConstrainedEnzymeIDs.append(enzymeID)


				else:
					enzymeEquivalentBoolID = self._generatedID_enzymeEquivBoolConstrained.format(enzymeID)

					self._solver.flowMaterialCoeffIs(
						reactionID,
						enzymeEquivalentBoolID,
						-1
						)


	def _initMass(self, externalExchangedMolecules, moleculeMasses):
		"""Create mass accumulation abstractions.

		Tracking the mass entering the system through metabolism is crucial for
		insuring closure. It can also be used to constrain the maximal rate of
		growth."""

		if moleculeMasses is not None:
			self._solver.flowMaterialCoeffIs(
				self._massOutName,
				self._massID,
				-1
				)

			for moleculeID in externalExchangedMolecules:
				exchangeFluxID = self._generatedID_externalExchange.format(moleculeID)

				try:
					moleculeMass = moleculeMasses[moleculeID]

				except KeyError:
					raise FBAError("You must provide masses for all molecules in externalExchangedMolecules")

				self._solver.flowMaterialCoeffIs(
					exchangeFluxID,
					self._massID,
					-moleculeMass # NOTE: negative because exchange fluxes point out
					)


	def _initMaintenance(self, maintenanceCost, maintenanceReaction):
		"""Create growth-associated maintenance abstractions.

		Two maintenance costs are typically associated with FBA; growth-
		associated maintenance is the energetic cost of increasing cell mass by
		a certain amount.  (Contrast non-growth-associated maintenance, which
		is a fixed energetic cost regardless of mass accumulation.)
		"""

		if (maintenanceCost is None) and (maintenanceReaction is None):
			return

		if (maintenanceCost is None) ^ (maintenanceReaction is None):
			raise FBAError("Must pass all or none of maintenanceCost, maintenanceReaction")


		# TODO: check that the mass flux stuff exists

		# computed mass output produces "GAM reactions"...
		reactionsNeededID = "GAM reactions" # TODO: move to class def

		self._solver.flowMaterialCoeffIs(
			self._massOutName,
			reactionsNeededID,
			maintenanceCost
			)

		# ... which are consumed in a seperate flux
		maintenanceReactionID = "Growth-associated maintenance" # TODO: move to class def

		self._solver.flowMaterialCoeffIs(
			maintenanceReactionID,
			reactionsNeededID,
			-1
			)

		for moleculeID, stoichCoeff in maintenanceReaction.viewitems():
			self._solver.flowMaterialCoeffIs(
				maintenanceReactionID,
				moleculeID,
				stoichCoeff
				)


	# Constraint setup

	# NOTE: I've used the nonspecific term "levels" since these could be
	# mole-concentration, number-concentration, or counts, depending on the
	# formulation of the problem.  All that matters is that initialization
	# parameters have consistent units.

	def externalMoleculeIDs(self):
		return self._externalMoleculeIDs


	def externalMoleculeLevelsIs(self, levels):
		levels_array = np.empty(len(self._externalMoleculeIDs))
		levels_array[:] = levels

		if (levels_array < 0).any():
			raise InvalidBoundaryError("Negative molecule levels not allowed")

		for moleculeID, level in izip(self._externalMoleculeIDs, levels_array):
			flowID = self._generatedID_externalExchange.format(moleculeID)

			self._solver.flowLowerBoundIs(
				flowID,
				-level
				)


	def internalMoleculeIDs(self):
		return self._internalMoleculeIDs


	def internalMoleculeLevelsIs(self, levels):
		levels_array = np.empty(len(self._internalMoleculeIDs))
		levels_array[:] = levels

		if (levels_array < 0).any():
			raise InvalidBoundaryError("Negative molecule levels not allowed")

		for moleculeID, level in izip(self._internalMoleculeIDs, levels_array):
			flowID = self._generatedID_internalExchange.format(moleculeID)

			self._solver.flowLowerBoundIs(
				flowID,
				-level
				)

			if self._forceInternalExchange:
				self._solver.flowUpperBoundIs(
					flowID,
					-level
					)


	def enzymeIDs(self):
		return self._enzymeIDs


	def enzymeLevelsIs(self, levels):
		if not hasattr(self, "_enzymeIDs"):
			return

		levels_array = np.empty(len(self._enzymeIDs))
		levels_array[:] = levels

		if (levels_array < 0).any():
			raise InvalidBoundaryError("Negative enzyme levels not allowed")

		for enzymeID, level in izip(self._enzymeIDs, levels_array):
			if enzymeID in self._rateConstrainedEnzymeIDs:
				# Rate-constrained
				flowID = self._generatedID_enzymeUsageRateConstrained.format(enzymeID)

				self._solver.flowUpperBoundIs(
					flowID,
					level
					)

			else:
				# Boolean-constrained (enzyme w/o an annotated rate)
				flowID = self._generatedID_enzymeUsageBoolConstrained.format(enzymeID)

				self._solver.flowUpperBoundIs(
					flowID,
					np.inf if level > 0 else 0
					)


	def reactionIDs(self):
		return np.array(self._reactionIDs)


	def maxReactionFluxIs(self, reactionID, maxFlux, raiseForReversible = True):
		if maxFlux < 0:
			raise InvalidBoundaryError("Maximum reaction flux must be at least 0")

		# if maxFlux < self._lowerBound[colIndex]:
		# 	raise InvalidBoundaryError("Maximum reaction flux must be greater than or equal to the minimum flux")

		reverseReactionID = self._generatedID_reverseReaction.format(reactionID)

		if raiseForReversible and reverseReactionID in self._reactionIDs:
			raise FBAError((
				"Setting the maximum reaction flux is ambiguous since " +
				"reaction {} has both a forward [{}] and reverse [{}] " +
				"component.  Call this method with argument " +
				"raiseForReversible = False if this is intended behavior."
				).format(reactionID, reactionID, reverseReactionID))

		self._solver.flowUpperBoundIs(
			reactionID,
			maxFlux
			)


	def minReactionFluxIs(self, reactionID, minFlux, raiseForReversible = True):
		if minFlux < 0:
			raise InvalidBoundaryError("Minimum reaction flux must be at least 0")

		# if minFlux > self._upperBound[colIndex]:
		# 	raise InvalidBoundaryError("Minimum reaction flux must be less than or equal to the maximum flux")

		reverseReactionID = self._generatedID_reverseReaction.format(reactionID)

		if raiseForReversible and reverseReactionID in self._reactionIDs:
			raise FBAError((
				"Setting the minimum reaction flux is ambiguous since " +
				"reaction {} has both a forward [{}] and reverse [{}] " +
				"component.  Call this method with argument " +
				"raiseForReversible = False if this is intended behavior."
				).format(reactionID, reactionID, reverseReactionID))

		self._solver.flowLowerBoundIs(
			reactionID,
			minFlux
			)

	# TODO: determine if this is needed

	# def objectiveIs(self, objective):
	# 	for moleculeID, coeff in objective.viewitems():
	# 		molecule_materialIndex = self._materialIndex(moleculeID)

	# 		pseudoFluxID = self._generatedID_moleculesToEquivalents.format(moleculeID)
	# 		colIndex = self._fluxIndex(pseudoFluxID)

	# 		self._A[molecule_materialIndex, colIndex] = -coeff


	def maxMassAccumulatedIs(self, maxAccumulation):
		self._solver.flowUpperBoundIs(
			self._massOutName,
			maxAccumulation
			)

	# Output

	def outputMoleculeIDs(self):
		return tuple(self._outputMoleculeIDs)


	def outputMoleculeLevelsChange(self):
		# This is essentially a dot product, need to profile to make sure this
		# isn't horribly slow

		change = np.zeros(len(self._outputMoleculeIDs))

		for i, outputMoleculeID in enumerate(self._outputMoleculeIDs):
			for reactionID, coeff in self._outputMoleculeCoeffs[i].viewitems():
				change[i] += self._solver.flowRates(reactionID) * coeff

		return -change


	# def externalExchangeFlux(self, moleculeID):
	# 	return -self._solver.flowRates(
	# 		self._generatedID_externalExchange.format(moleculeID)
	# 		)


	def externalExchangeFluxes(self):
		return self._solver.flowRates(self._externalExchangeIDs)


	# def internalExchangeFlux(self, moleculeID):
	# 	# TODO
	# 	return -self._solutionFluxes[
	# 		self._fluxIndex(self._generatedID_internalExchange.format(moleculeID))
	# 		]


	def reactionFlux(self, reactionID):
		return self._solver.flowRates(reactionID)


	def reactionFluxes(self):
		return self._solver.flowRates(self._reactionIDs)


	def objectiveReactionFlux(self): # TODO: rename to biomassReactionFlux
		# catch exceptions
		return self._solver.flowRates(self._standardObjectiveReactionName)


	# def objectiveValue(self):
	# 	return self._rawSolution["primal objective"]


	# def enzymeUsage(self):
	# 	return self._solutionFluxes[self._enzymeUsageRateConstrainedIndexes]


	# def massAccumulated(self):
	# 	return self._solutionFluxes[self._fluxIndex(self._massOutName)]
