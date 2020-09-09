#!/usr/bin/env python

"""
EnzymeKinetics

Takes in enzyme kinetics data on initialization, and returns dicts of rate estimates when passed
metabolite and enzyme concentrations at runtime.
"""
from __future__ import absolute_import, division, print_function

import numpy as np

from wholecell.utils import units
from Equation import Expression
import six

class enzymeKineticsError(Exception):
	pass

COUNTS_UNITS = units.umol
TIME_UNITS = units.s
VOLUME_UNITS = units.L

class EnzymeKinetics(object):
	"""
	EnzymeKinetics

	Returns rate estimates from kinetic equation information stored in reactionRateInfo.
	"""

	def __init__(self, reactionRateInfo, kcatsOnly=False, useCustoms=True, moreThanKcat=False):

		# Set default reaction rate limit, to which reactions are set absent other information
		self.defaultRate = (COUNTS_UNITS / TIME_UNITS / VOLUME_UNITS) * np.inf

		# Load rate functions from enzymeKinetics.tsv flat file
		self.reactionRateInfo = reactionRateInfo

		## Filter the reactions as specified
		# Exclude any rows with more than a kcat
		if kcatsOnly:
			reactionRateInfoNew = {}
			for constraintID, reactionInfo in six.viewitems(self.reactionRateInfo):
				# Kcat-only reactions will have no kMs, kIs, or custom equations
				if len(reactionInfo["kM"]) or len(reactionInfo["kI"]) or reactionInfo["customRateEquation"]:
					continue
				reactionRateInfoNew[constraintID] = reactionInfo
			self.reactionRateInfo = reactionRateInfoNew

		# Exclude any custom equation rows
		if not useCustoms:
			reactionRateInfoNew = {}
			for constraintID, reactionInfo in six.viewitems(self.reactionRateInfo):
				if reactionInfo["customRateEquation"] is None:
					reactionRateInfoNew[constraintID] = reactionInfo
			self.reactionRateInfo = reactionRateInfoNew

		# Throw out any kcat-only reactions
		if moreThanKcat:
			reactionRateInfoNew = {}
			for constraintID, reactionInfo in six.viewitems(self.reactionRateInfo):
				if len(reactionInfo["kM"]) or len(reactionInfo["kI"]) or reactionInfo["customRateEquation"]:
					reactionRateInfoNew[constraintID] = reactionInfo
			self.reactionRateInfo = reactionRateInfoNew

		self.allConstraintIDs = list(self.reactionRateInfo.keys())

		self.allReactionIDs = [x["reactionID"] for x in self.reactionRateInfo.values()]

		self.inputsChecked = False

	def checkKnownSubstratesAndEnzymes(self, metaboliteSMatrixNamesNoCompartment, metaboliteNamesWithConcentrations, enzymeNames, removeUnknowns=False):
		knownConstraints = {}
		unusableConstraints = {}
		unknownSubstrates = set()
		unknownEnzymes = set()
		unknownCustomVars = set()


		for constraintID, reactionInfo in six.viewitems(self.reactionRateInfo):
			keepReaction = True
			reactionType = reactionInfo["rateEquationType"]
			if reactionType == "standard":

				# Check if the substrates are known
				if len(reactionInfo["kM"]) > 0:
					for substrateID in reactionInfo["substrateIDs"]:
						if substrateID not in metaboliteNamesWithConcentrations:
							unknownSubstrates.add(substrateID)
							unusableConstraints[constraintID] = reactionInfo
							keepReaction = False
				else:
					for substrateID in reactionInfo["substrateIDs"]:
						if substrateID[:-3] not in metaboliteSMatrixNamesNoCompartment:
							unknownSubstrates.add(substrateID)
							unusableConstraints[constraintID] = reactionInfo
							keepReaction = False


				# Check if the enzymes are known
				for enzymeID in reactionInfo["enzymeIDs"]:
					if enzymeID not in enzymeNames:
						unknownEnzymes.add(enzymeID)
						unusableConstraints[constraintID] = reactionInfo
						keepReaction = False


			elif reactionType == "custom":

				for variable in reactionInfo["customParameterVariables"].values():
					if variable not in metaboliteNamesWithConcentrations:
						notSubstrate = True
					if variable not in enzymeNames:
						notEnzyme = True

					if notSubstrate and notEnzyme:
						unknownCustomVars.add(variable)
						unusableConstraints[constraintID] = reactionInfo
						keepReaction = False
			else:
				# Reaction type is unknown
				raise Exception("Reaction type '%s' is unknown. Must be either 'standard' or 'custom'." % (reactionType,))

			# Keep the reaction only if both substrates and enzymes are known
			if keepReaction:
				knownConstraints[constraintID] = reactionInfo

		if removeUnknowns:
			self.reactionRateInfo = knownConstraints
			self.inputsChecked = True

		unknownVals = {"unknownSubstrates":unknownSubstrates, "unknownEnzymes":unknownEnzymes, "unknownCustomVars":unknownCustomVars}

		return knownConstraints, unusableConstraints, unknownVals


	def reactionRate(self, reactionInfo, metaboliteConcentrationsDict, enzymeConcentrationsDict):

		if reactionInfo["rateEquationType"] == "standard":
			return self.reactionStandard(reactionInfo, metaboliteConcentrationsDict, enzymeConcentrationsDict)
		elif reactionInfo["rateEquationType"] == "custom":
			return self.reactionCustom(reactionInfo, metaboliteConcentrationsDict, enzymeConcentrationsDict)
		else:
			raise NameError("rateEquationType %s not recognized! Must be either 'standard' or 'custom'." % (reactionInfo["rateEquationType"]))

	def reactionStandard(self, reactionInfo, metaboliteConcentrationsDict, enzymeConcentrationsDict):
		# Find the enzymes needed for this rate
		enzymeConc = enzymeConcentrationsDict[reactionInfo["enzymeIDs"][0]].asNumber(COUNTS_UNITS/VOLUME_UNITS)
		kMs = reactionInfo["kM"]
		kIs = reactionInfo["kI"]
		substrateIDs = reactionInfo["substrateIDs"]

		if len(kMs) + len(kIs) > len(substrateIDs):
			raise enzymeKineticsError("The number of saturation constants (kMs and kIs) must not be greater than the number of substrates for a standard reaction. For constraint {}, there are {} kMs, {} kIs, and {} substrates. ReactionInfo: {}".format(reactionInfo["constraintID"], len(kMs), len(kIs), len(substrateIDs), reactionInfo))

		rate = np.amax(reactionInfo["kcat"]) * enzymeConc
		if "kcatAdjusted" in reactionInfo:
			rate = np.amax(reactionInfo["kcatAdjusted"]) * enzymeConc

		idx = 0
		for kM in reactionInfo["kM"]:
			substrateConc = metaboliteConcentrationsDict[substrateIDs[idx]].asNumber(COUNTS_UNITS/VOLUME_UNITS)
			rate *= (substrateConc / (float(substrateConc) + kM))
			idx+=1
		for kI in reactionInfo["kI"]:
			substrateConc = metaboliteConcentrationsDict[substrateIDs[idx]].asNumber(COUNTS_UNITS/VOLUME_UNITS)
			rate *= 1.0 / (1.0 + (substrateConc / kI))
			idx+=1

		return (COUNTS_UNITS / TIME_UNITS / VOLUME_UNITS) * rate

	def reactionCustom(self, reactionInfo, metaboliteConcentrationsDict, enzymeConcentrationsDict):
		enzymeConcArray = [enzymeConcentrationsDict[reactionInfo["enzymeIDs"][0]]].asNumber(COUNTS_UNITS/VOLUME_UNITS)

		customParameters = reactionInfo["customParameters"]
		customParameterVariables = reactionInfo["customParameterVariables"]
		customParameterConstants = reactionInfo["customParameterConstantValues"]
		equationString = reactionInfo["customRateEquation"]
		parameterDefinitionArray = reactionInfo["customParameters"]

		parametersArray = customParameterConstants
		for customParameter in customParameters[len(customParameterConstants):]:
			variable = customParameterVariables[customParameter]
			if variable in enzymeConcentrationsDict:
				parametersArray.append(enzymeConcentrationsDict[variable].asNumber(COUNTS_UNITS/VOLUME_UNITS))
			elif variable in metaboliteConcentrationsDict:
				parametersArray.append(metaboliteConcentrationsDict[variable].asNumber(COUNTS_UNITS/VOLUME_UNITS))

		assert (len(parametersArray) == len(parameterDefinitionArray))

		customRateLaw = Expression(equationString, parameterDefinitionArray)

		return (COUNTS_UNITS / TIME_UNITS / VOLUME_UNITS) * customRateLaw(*parametersArray)


	def allConstraintsDict(self, metaboliteConcentrationsDict, enzymeConcentrationsDict):
		constraintsDict = {}
		for constraintID, reactionInfo in six.viewitems(self.reactionRateInfo):
			constraintsDict[constraintID] = self.reactionRate(reactionInfo, metaboliteConcentrationsDict, enzymeConcentrationsDict)

		return constraintsDict

	def allReactionsDict(self, metaboliteConcentrationsDict, enzymeConcentrationsDict):
		"""
		Create a dict of dicts from reactionID to constraintIDs for that reaction, to rates for each constraintID.
		"""
		reactionsDict = {}
		for constraintID, reactionInfo in six.viewitems(self.reactionRateInfo):
			reactionID = reactionInfo["reactionID"]
			if reactionID not in reactionsDict:
				reactionsDict[reactionID] = {}
			reactionsDict[reactionID][constraintID] = self.reactionRate(reactionInfo, metaboliteConcentrationsDict, enzymeConcentrationsDict)

		return reactionsDict

	def ratesView(self, reactionIDs, reactionsToConstraintsDict, metaboliteConcentrationsDict, enzymeConcentrationsDict, raiseIfNotFound=False):
		"""
		Returns an array of rates with units.
			Order taken from reactionIDs, rates to estimate come from reactionsToConstraintsDict.
			When a rate is not found, raises exception if raiseIfNotFound, else returns default rate.
			A reaction must be in reactionsToConstraintsDict, or it will get the default value.
		"""

		unknownConstraints = set()
		unknownReactions = set()

		# Build an estimate for rates of constraints passed in
		rates = self.defaultRate * np.ones(len(reactionIDs))
		for idx, reactionID in enumerate(reactionIDs):
			if reactionID in reactionsToConstraintsDict:
				constraintID = reactionsToConstraintsDict[reactionID]["constraintID"]
				coefficient = reactionsToConstraintsDict[reactionID]["coefficient"] if "coefficient" in reactionsToConstraintsDict[reactionID] else 1
				if constraintID in self.reactionRateInfo:
					rates[idx] = coefficient * self.reactionRate(self.reactionRateInfo[constraintID], metaboliteConcentrationsDict, enzymeConcentrationsDict)
				else:
					if raiseIfNotFound:
						unknownConstraints.add(constraintID)
			else:
				if raiseIfNotFound:
					unknownReactions.add(reactionID)

		if len(unknownConstraints) > 0:
			raise Exception("No rate estimate found for the following {} constraintIDs {}.".format(len(unknownConstraints), unknownConstraints))

		if len(unknownReactions) > 0:
			raise Exception("No rate estimate found for the following {} reactionIDs {}.".format(len(unknownReactions), unknownReactions))

		return rates

	def ratesViewConstraints(self, constraintIDs, metaboliteConcentrationsDict, enzymeConcentrationsDict, raiseIfNotFound=False):
		"""
		Returns an array of rates with units, in the same order as the constraintIDs iterable.
		"""
		# Check if all needed metabolite and enzyme concentrations are given
		if not self.inputsChecked:
			knownConstraints, unusableConstraints, unknownVals = self.checkKnownSubstratesAndEnzymes(metaboliteConcentrationsDict, enzymeConcentrationsDict, enzymeNames={}, removeUnknowns=False)
			if len(unusableConstraints) > 0:
				raise Exception("Unable to compute kinetic rate for these reactions: {}\n. Missing values for: {}".format(list(unusableConstraints.keys()), unknownVals))

		unknownConstraints = set()

		# Build an estimate for rates of constraints passed in
		rates = self.defaultRate * np.ones(len(constraintIDs))
		for idx, constraintID in enumerate(constraintIDs):
			if constraintID in self.reactionRateInfo:
				rates[idx] = self.reactionRate(self.reactionRateInfo[constraintID], metaboliteConcentrationsDict, enzymeConcentrationsDict)
			else:
				if raiseIfNotFound:
					unknownConstraints.add(constraintID)

		if len(unknownConstraints) > 0:
			raise Exception("No rate estimate found for constraintIDs {}.".format(unknownConstraints))

		return rates


	def ratesViewReactions(self, reactionIDs, metaboliteConcentrationsDict, enzymeConcentrationsDict, sortFunction, raiseIfNotFound=False):
		"""
		Returns an array of rates with units, in the same order as the reactionIDs.
		Uses sortFunction to decide which constraint to use if a reaction has multiple choices.
		"""
		# Check if all needed metabolite and enzyme concentrations are given
		if not self.inputsChecked:
			knownConstraints, unusableConstraints, unknownVals = self.checkKnownSubstratesAndEnzymes(metaboliteConcentrationsDict, enzymeConcentrationsDict, enzymeNames={}, removeUnknowns=False)
			if len(unusableConstraints) > 0:
				raise Exception("Unable to compute kinetic rate for these reactions: {}\n. Missing values for: {}".format(list(unusableConstraints.keys()), unknownVals))

		unknownReactions = set()

		rates = self.defaultRate * np.ones(len(reactionIDs))

		constraintsDict = self.allConstraintsDict(metaboliteConcentrationsDict, enzymeConcentrationsDict)
		for idx, reactionID in enumerate(reactionIDs):
			if reactionID in constraintsDict:
				rates[idx] = sortFunction(constraintsDict[reactionID])
			else:
				if raiseIfNotFound:
					unknownReactions.add(reactionID)

		if len(unknownReactions) > 0:
			raise Exception("No rate estimate found for reactionIDs {}.".format(unknownReactions))


		return rates

	def ratesViewLowest(self, reactionIDs, metaboliteConcentrationsDict, enzymeConcentrationsDict, raiseIfNotFound=False):
		"""
		Returns an array of rates with units, in the same order as the reactionIDs, chosing the minimum estimate when a reactionID has multiple constraintIDs.
		"""
		return self.ratesViewReactions(reactionIDs, metaboliteConcentrationsDict, enzymeConcentrationsDict, np.amin, raiseIfNotFound)


	def ratesViewHighest(self, reactionIDs, metaboliteConcentrationsDict, enzymeConcentrationsDict, raiseIfNotFound=False):
		"""
		Returns an array of rates with units, in the same order as the reactionIDs, chosing the minimum estimate when a reactionID has multiple constraintIDs.
		"""
		return self.ratesViewReactions(reactionIDs, metaboliteConcentrationsDict, enzymeConcentrationsDict, np.amax, raiseIfNotFound)