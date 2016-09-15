"""
SimulationData for metabolism process

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 03/06/2015
"""

from __future__ import division

from itertools import chain

from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray
import numpy as np
import collections
import warnings

ILE_LEU_CONCENTRATION = 3.0e-4 # mmol/L
ILE_FRACTION = 0.360 # the fraction of iso/leucine that is isoleucine; computed from our monomer data
ECOLI_PH = 7.2

PPI_CONCENTRATION = 0.5e-3 # M, multiple sources

EXCHANGE_UNITS = units.mmol / units.g / units.h

# If true, enzyme kinetics entries which reference unknown reactions are ignored
# If false, raises an exception in such a case
raiseForUnknownRxns = False
warnForUnknownRxns = False
raiseForTruncatedRxns = False
warnForTruncatedRxns = False
warnForMultiTruncatedRxns = False

reverseReactionString = "{} (reverse)"


class Metabolism(object):
	""" Metabolism """

	def __init__(self, raw_data, sim_data):
		self._buildBiomass(raw_data, sim_data)
		self._buildMetabolism(raw_data, sim_data)


	def _buildBiomass(self, raw_data, sim_data):
		wildtypeIDs = set(entry["molecule id"] for entry in raw_data.biomass)
		# TODO: unjank this

		# Load the biomass function flat file as a dict
		self.biomassFunction = {entry['molecule id']:entry['coefficient'] for entry in raw_data.biomass}

		self.previousBiomassMeans = {entry['molecule id']:entry['mean flux'] for entry in raw_data.previousBiomassFluxes}
		self.previousBiomassLog10Means = {entry['molecule id']:entry['mean log10 flux'] for entry in raw_data.previousBiomassFluxes}
		self.previousBiomassStds = {entry['molecule id']:entry['standard deviation'] for entry in raw_data.previousBiomassFluxes}

		# Create vector of metabolite pools (concentrations)

		# Since the data only covers certain metabolites, we need to rationally
		# expand the dataset to include the other molecules in the biomass
		# function.

		# First, load in metabolites that do have concentrations, then assign
		# compartments according to those given in the biomass objective.  Or,
		# if there is no compartment, assign it to the cytoplasm.

		metaboliteIDs = []
		metaboliteConcentrations = []
		metaboliteConcentrationData = dict((m["Metabolite"], m["Concentration"].asNumber(units.mol / units.L)) for m in raw_data.metaboliteConcentrations)

		wildtypeIDtoCompartment = {
			wildtypeID[:-3] : wildtypeID[-3:]
			for wildtypeID in wildtypeIDs
			} # this assumes biomass reaction components only exist in a single compartment

		for metaboliteID, concentration in metaboliteConcentrationData.viewitems():
			if metaboliteID in wildtypeIDtoCompartment:
				metaboliteIDs.append(
					metaboliteID + wildtypeIDtoCompartment[metaboliteID]
					)

			else:
				metaboliteIDs.append(
					metaboliteID + "[c]"
					)

			metaboliteConcentrations.append(concentration)

		# ILE/LEU: split reported concentration according to their relative abundances
		# TODO: more thorough estimate of abundance or some external data point (neidhardt?)

		ileRelative = ILE_FRACTION
		leuRelative = 1 - ileRelative

		metaboliteIDs.append("ILE[c]")
		metaboliteConcentrations.append(ileRelative * ILE_LEU_CONCENTRATION)

		metaboliteIDs.append("LEU[c]")
		metaboliteConcentrations.append(leuRelative * ILE_LEU_CONCENTRATION)

		# CYS/SEC/GLY: fit a relative abundance:concentration line (L1 norm)
		# with other amino acids and solve for these

		aaConcentrations = []
		# aaAbundancesWithConcentrations = []

		for aaIndex, aaID in enumerate(sim_data.amino_acid_1_to_3_ordered.values()):
			if aaID in metaboliteIDs:
				metIndex = metaboliteIDs.index(aaID)
				aaConcentrations.append(metaboliteConcentrations[metIndex])

		# TODO: implement L1-norm minimization

		# for now: just choosing and assigning the smallest value
		aaSmallestConc = min(aaConcentrations)

		# HACK: min conc. doesn't work here
		metaboliteIDs.append("GLY[c]")
		metaboliteConcentrations.append(
			metaboliteConcentrations[metaboliteIDs.index("L-ALPHA-ALANINE[c]")]
			)

		metaboliteIDs.append("CYS[c]")
		metaboliteConcentrations.append(aaSmallestConc)

		metaboliteIDs.append("L-SELENOCYSTEINE[c]")
		metaboliteConcentrations.append(aaSmallestConc)

		# DGTP: set to smallest of all other DNTP concentrations

		dntpConcentrations = []
		# dntpAbundancesWithConcentrations = []

		for dntpIndex, dntpID in enumerate(sim_data.moleculeGroups.dNtpIds):
			if dntpID in metaboliteIDs:
				metIndex = metaboliteIDs.index(dntpID)
				dntpConcentrations.append(metaboliteConcentrations[metIndex])

		dntpSmallestConc = min(dntpConcentrations)

		metaboliteIDs.append("DGTP[c]")
		metaboliteConcentrations.append(dntpSmallestConc)

		# H: reported pH

		hydrogenConcentration = 10**(-ECOLI_PH)

		metaboliteIDs.append("PROTON[c]")
		metaboliteConcentrations.append(hydrogenConcentration)

		# PPI: multiple sources report 0.5 mM

		# NOTE: Nick says that the physiological levels of PPI are very low - investigate this

		metaboliteIDs.append("PPI[c]")
		metaboliteConcentrations.append(PPI_CONCENTRATION)

		metaboliteIDs.append("Pi[c]") # TODO: find a real number
		metaboliteConcentrations.append(PPI_CONCENTRATION)

		# NOTE: this assertion is thrown since there are many biomass things not in the (metabolic) model
		# unaccounted = set(wildtypeIDs) - set(metaboliteIDs)
		# assert len(unaccounted) == 0

		# Add byproducts with no annotated concentration to force recycling

		metaboliteIDs.append("UMP[c]")
		metaboliteConcentrations.append(2.40e-5)

		# Other quantities to consider:
		# - (d)NTP byproducts not currently included

		for key, value in sim_data.mass.getBiomassAsConcentrations(sim_data.doubling_time).iteritems():
			metaboliteIDs.append(key)
			metaboliteConcentrations.append(value.asNumber(units.mol / units.L))

		self.nutrientData = sim_data.nutrientData
		self.concentrationUpdates = ConcentrationUpdates(dict(zip(
			metaboliteIDs,
			(units.mol / units.L) * np.array(metaboliteConcentrations)
			)),
			raw_data.equilibriumReactions,
			self.nutrientData,
		)
		self.concDict = self.concentrationUpdates.concentrationsBasedOnNutrients("minimal")
		self.nutrientsToInternalConc = {}
		self.nutrientsToInternalConc["minimal"] = self.concDict.copy()


	def _buildMetabolism(self, raw_data, sim_data):
		# Build the matrices/vectors for metabolism (FBA)

		# These may be modified/extended later, but should provide the basic
		# data structures

		reactionStoich = {}
		reversibleReactions = []
		reactionEnzymes = {}
		reactionRates = {}

		enzymeExceptions = set()

		validEnzymeIDs = set([])
		validProteinIDs = ['{}[{}]'.format(x['id'],location) for x in raw_data.proteins for location in x['location']]
		validProteinComplexIDs = ['{}[{}]'.format(x['id'],location) for x in raw_data.proteinComplexes for location in x['location']]
		validEnzymeIDs.update(validProteinIDs)
		validEnzymeIDs.update(validProteinComplexIDs)
		validEnzymeCompartments = collections.defaultdict(set)

		self.default_kcat = raw_data.parameters["carbonicAnhydraseKcat"]

		for enzymeID in validEnzymeIDs:
			enzyme = enzymeID[:enzymeID.index("[")]
			location = enzymeID[enzymeID.index("[")+1:enzymeID.index("[")+2]

			validEnzymeCompartments[enzyme].add(location)

		# Enzymes which should not be used for enzyme-reaction pairs
		for rxnEnzymePair in raw_data.unconstrainedReactionEnzymes:
			enzymeExceptions.add(rxnEnzymePair["enzymeID"])

		for reaction in raw_data.reactions:
			reactionID = reaction["reaction id"]
			stoich = reaction["stoichiometry"]
			reversible = reaction["is reversible"]
			enzyme_list = self.addEnzymeCompartmentTags(reaction["catalyzed by"], validEnzymeCompartments)

			if len(stoich) <= 1:
				raise Exception("Invalid biochemical reaction: {}, {}".format(reactionID, stoich))

			reactionStoich[reactionID] = stoich

			# Remove enzyme-reaction links for any enzyme in enzymeExceptions
			for enzymeID in enzyme_list:
				if enzymeID in enzymeExceptions:
					enzyme_list.remove(enzymeID)
			enzymeKcatLink = {enzymeID:self.default_kcat.asNumber(1 / units.s) for enzymeIDs in enzyme_list}
			reactionEnzymes[reactionID] = enzymeKcatLink

			# Add the reverse reaction
			if reversible:
				reverseReactionID = reverseReactionString.format(reactionID)
				reactionStoich[reverseReactionID] = {
					moleculeID:-stoichCoeff
					for moleculeID, stoichCoeff in reactionStoich[reactionID].viewitems()
					}

				reactionEnzymes[reverseReactionID] = enzymeKcatLink
				reversibleReactions.append(reactionID)

		reactionRateInfo = {}
		constraintIDs = []
		constraintToReactionDict = {}

		directionAmbiguousRxns = set()
		directionInferedReactions = set()
		nonCannonicalRxns = set()
		unknownRxns = set()
		truncatedRxns = {}
		multipleTruncatedOptionRxns = collections.defaultdict(set)

		# Enzyme kinetics data
		for idx, reaction in enumerate(raw_data.enzymeKinetics):
			reactionID = reaction["reactionID"]
			# Ensure all reactions in enzymeKinetics refer to tha actual corresponding reaction in reactions.tsv, rather than a non-canonical alternative substrate
			if reactionID not in reactionStoich:
				truncated = False
				for reactionName in reactionStoich:
					if reactionName.startswith(reactionID):
						truncated = True
						if reactionName in truncatedRxns:
							if truncatedRxns[reactionID] != reactionName:
								multipleTruncatedOptionRxns[reactionID].add(truncatedRxns[reactionID])
								del truncatedRxns[reactionID]
								multipleTruncatedOptionRxns[reactionID].add(reactionName)
						elif reactionName in multipleTruncatedOptionRxns:
							multipleTruncatedOptionRxns[reactionID].add(reactionName)
						else:
							truncatedRxns[reactionID] = reactionName
				if not truncated:
					unknownRxns.add(reaction["reactionID"])
				continue

		for idx, reaction in enumerate(raw_data.enzymeKinetics):
			reactionID = reaction["reactionID"]

			# Find the temp-adjusted kcat value
			newKcat = self.temperatureAdjustedKcat(reaction)

			# Add compartment tags to enzymes
			reaction["enzymeIDs"] = self.addEnzymeCompartmentTags(reaction["enzymeIDs"], validEnzymeCompartments)

			# If the substrates don't already have a compartment tag, add [c] (cytosol) as a default
			reaction["substrateIDs"] = [x + '[c]' if x[-3:-2] != '[' else x for x in reaction["substrateIDs"]]

			if reaction["rateEquationType"] == 'custom':
				# If the enzymes and substrates mentioned in the
				# customParameterVariables dict lack a compartment tag, add [c]
				# (cytosol) as a default
				parametersDict = reaction["customParameterVariables"]
				for key in parametersDict:
					value = parametersDict[key]
					if value[-3:-2] != '[':
						parametersDict[key] = value + '[c]'
				reaction["customParameterVariables"] = parametersDict

			# Ensure all reactions in enzymeKinetics refer to tha actual corresponding reaction in reactions.tsv, rather than a non-canonical alternative substrate
			if reactionID in reactionStoich:
				thisRxnStoichiometry = reactionStoich[reactionID]
			else:
				if reactionID in truncatedRxns:
					thisRxnStoichiometry = reactionStoich[truncatedRxns[reactionID]]
				elif reactionID in multipleTruncatedOptionRxns:
					reactionOptions = set(multipleTruncatedOptionRxns[reactionID])

					# Infer on reaction substrates
					kineticSubstrates = reaction["substrateIDs"]
					threshold = 0
					candidates = set()
					ambiguous = False
					for reactionOption in reactionOptions:
						substrates = set(reactionStoich[reactionOption].keys())
						amountOfOverlap = len(substrates.intersection(kineticSubstrates))
						if amountOfOverlap > threshold:
							threshold = amountOfOverlap
							candidates = set([reactionOption])
							ambiguous = False
						elif amountOfOverlap == threshold:
							if threshold > 0:
								ambiguous = True
							candidates.add(reactionOption)
					if threshold == 0:
						# Take none, add this reaction to the unknown reactions list, remove it from multipleTruncatedOptionRxns
						unknownRxns.add(reactionID)
						del multipleTruncatedOptionRxns[reactionID]
					if ambiguous:
						# Keep all highest-scoring reactions in multipleTruncatedOptionRxns, remove any reactions which score less than threshold
						multipleTruncatedOptionRxns[reactionID] =  candidates
					else:
						# Take the highest-scoring reaction
						assert len(candidates) == 1
						thisRxnStoichiometry = reactionStoich[candidates.pop()]

						# Remove this reaction from multipleTruncatedOptionRxns
						del multipleTruncatedOptionRxns[reactionID]

				elif reactionID in unknownRxns:
					continue
				else:
					raise Exception("Something went wrong with reaction {}".format(reactionID))

			substrateIDs = reaction["substrateIDs"]
			if reaction["rateEquationType"] == "standard":
				if len(reaction["kI"]) > 0:
					for substrate in substrateIDs[:-len(reaction["kI"])]:
						if substrate not in thisRxnStoichiometry.keys():
							nonCannonicalRxns.add(reaction["constraintID"])
			elif reaction["rateEquationType"] == "custom":
				pass
			else:
				raise Exception("rateEquationType {} not understood in reaction {} on enzymeKinetics line {}".format(reaction["reactionID"], reaction["reactionID"], idx))

			# Check if this constraint is for a reverse reaction
			if reactionID in reversibleReactions:
				if reaction["direction"] == "forward":
					pass
				elif reaction["direction"] == "reverse":
					reaction["reactionID"] = reverseReactionString.format(reactionID)
					reaction["constraintID"] = reverseReactionString.format(reaction["constraintID"])
				else:
					# Infer directionality from substrates
					directionInferedReactions.add(reaction["reactionID"])
					if reaction["rateEquationType"] == "standard":
						reverseCounter = 0.
						allCounter = 0.
						# How many substrates are products, how many are reactants?
						for substrate in reaction["substrateIDs"]:
							if substrate in thisRxnStoichiometry.keys():
								if thisRxnStoichiometry[substrate] > 0:
									reverseCounter += 1.
							allCounter += 1.
						# If more are products than reactants, treat this as a reverse reaction constraint
						# Record any ambiguous reactions and throw an exception once all are gathered
						if allCounter < 1.0:
							directionAmbiguousRxns.add(reaction["reactionID"])
							continue

						if reverseCounter == allCounter:
							reaction["reactionID"] = reverseReactionString.format(reactionID)
							reaction["constraintID"] = reverseReactionString.format(reaction["constraintID"])
						elif reverseCounter > 0:
							if len(reaction["kI"]) == reverseCounter:
								pass
							else:
								directionAmbiguousRxns.add(reaction["reactionID"])
					elif reaction["rateEquationType"] == "custom":
						raise Exception("Custom equations for reversible reactions must specify direction. Reaction {} on enzymeKinetics line {} does not.".format(reactionID, idx))
					else:
						raise Exception("Reaction {} with rateEquationType ({}) not recognized on enzymeKinetics line {}.".format(reactionID, reaction["rateEquationType"], idx))

			constraintID = reaction["constraintID"]
			constraintIDs.append(constraintID)
			constraintToReactionDict[constraintID] = reactionID

			reactionRateInfo[constraintID] = reaction

		if len(unknownRxns) > 0:
			message = "The following {} enzyme kinetics reactions appear to be for reactions which don't exist in the model - they should be corrected or removed. {}".format(len(unknownRxns), unknownRxns)
			if raiseForUnknownRxns:
				raise Exception(message)
			elif warnForUnknownRxns:
				warnings.warn(message)

		if len(truncatedRxns) > 0:
			message = "The following {} enzyme kinetics reaction names are truncated versions of reactions in the model {}".format(len(truncatedRxns), truncatedRxns)
			if warnForTruncatedRxns:
				warnings.warn(message)
			elif raiseForTruncatedRxns:
				raise Exception(message)

		if len(multipleTruncatedOptionRxns) > 0:
			message = "The following {} enzyme kinetics reaction names are truncated versions of more than one reaction in the model {}".format(len(multipleTruncatedOptionRxns), multipleTruncatedOptionRxns)
			if warnForMultiTruncatedRxns:
				warnings.warn(message)
			elif raiseForTruncatedRxns:
				raise Exception(message)

		if len(directionAmbiguousRxns) > 0:
			raise Exception("The following enzyme kinetics entries have ambiguous direction. Split them into multiple lines in the flat file to increase clarity. {}".format(directionAmbiguousRxns))

		if len(nonCannonicalRxns) > 0:
			raise Exception("The following {} enzyme kinetics entries reference substrates which don't appear in their corresponding reaction, and aren't paired with an inhibitory constant (kI). They should be corrected or removed. {}".format(len(nonCannonicalRxns), nonCannonicalRxns))

		self.reactionEnzymes = self.buildEnzymeReactionKcatLinks(reactionRateInfo, reactionEnzymes)

		self.reactionStoich = reactionStoich
		self.nutrientsTimeSeries = sim_data.nutrientsTimeSeries
		self.maintenanceReaction = {"ATP[c]": -1, "WATER[c]": -1, "ADP[c]": +1, "Pi[c]": +1, "PROTON[c]": +1,}
		self.reversibleReactions = reversibleReactions
		self.directionInferedReactions = sorted(list(directionInferedReactions))
		self.reactionRateInfo = reactionRateInfo
		self.enzymeNames = list(validEnzymeIDs)
		self.constraintIDs = constraintIDs
		self.constraintToReactionDict = constraintToReactionDict

	def exchangeConstraints(self, exchangeIDs, coefficient, targetUnits, nutrientsTimeSeriesLabel, time, concModificationsBasedOnCondition = None, preview = False):
		newObjective = None
		while len(self.nutrientsTimeSeries[nutrientsTimeSeriesLabel]) and time > self.nutrientsTimeSeries[nutrientsTimeSeriesLabel][0][0]:
			if preview:
				_, nutrients = self.nutrientsTimeSeries[nutrientsTimeSeriesLabel][0]
			else:
				_, nutrients = self.nutrientsTimeSeries[nutrientsTimeSeriesLabel].popleft()
			self._unconstrainedExchangeMolecules = self.nutrientData["importUnconstrainedExchangeMolecules"][nutrients]
			self._constrainedExchangeMolecules = self.nutrientData["importConstrainedExchangeMolecules"][nutrients]
			concDict = self.concentrationUpdates.concentrationsBasedOnNutrients(nutrients, self.nutrientsToInternalConc)
			if concModificationsBasedOnCondition is not None:
				concDict.update(concModificationsBasedOnCondition)
			newObjective = dict((key, concDict[key].asNumber(targetUnits)) for key in concDict)
			if preview:
				break

		externalMoleculeLevels = np.zeros(len(exchangeIDs), np.float64)

		for index, moleculeID in enumerate(exchangeIDs):
			if moleculeID in self._unconstrainedExchangeMolecules:
				externalMoleculeLevels[index] = np.inf

			elif moleculeID in self._constrainedExchangeMolecules.viewkeys():
				externalMoleculeLevels[index] = (
					self._constrainedExchangeMolecules[moleculeID] * coefficient
					).asNumber(targetUnits)
			else:
				externalMoleculeLevels[index] = 0.

		return externalMoleculeLevels, newObjective

	def addEnzymeCompartmentTags(self, enzymeIDs, validEnzymeCompartmentsDict):
		"""
		If the enzymes don't already have a compartment tag, add one from the valid compartment list or [c] (cytosol) as a default
		"""
		new_reaction_enzymes = []
		for reactionEnzyme in enzymeIDs:
			if reactionEnzyme[-3:-2] !='[':
				if len(validEnzymeCompartmentsDict[reactionEnzyme]) > 0:
					new_reaction_enzymes.append(reactionEnzyme +'['+str(validEnzymeCompartmentsDict[reactionEnzyme].pop())+']')
				else:
					new_reaction_enzymes.append(reactionEnzyme + '[c]')
			else:
				new_reaction_enzymes.append(reactionEnzyme)

		return new_reaction_enzymes

	def enzymeReactionMatrix(self, reactionIDs, enzymeNames, reactionEnzymesDict):
		"""
		Builds a (num reactions) by (num enzymes) matrix which maps enzyme concentrations to overall reaction rate.
		reactionEnzymesDict is a dict from reactionID:{dict of enzymes catalyzing this reaction:their associated kcat}
		"""
		assert sorted(reactionIDs) == sorted(reactionEnzymesDict.keys())

		enzymeNames = list(enzymeNames)

		enzymeReactionMatrix = np.zeros((len(reactionIDs),len(enzymeNames)))
		for rxnIdx, reactionID in enumerate(reactionIDs):
			if reactionID in reactionEnzymesDict:
				for enzymeName, kcat in reactionEnzymesDict[reactionID].iteritems():
					if enzymeName in enzymeNames:
						enzymeIdx = enzymeNames.index(enzymeName)
						enzymeReactionMatrix[rxnIdx, enzymeIdx] = kcat
		return enzymeReactionMatrix

	def buildEnzymeReactionKcatLinks(self, reactionRateInfo, reactionEnzymesDict):
		for constraintID, reactionInfo in reactionRateInfo.iteritems():
			if reactionInfo["rateEquationType"] == "custom":
				continue

			reactionID = reactionInfo["reactionID"]
			enzymeIDs = reactionInfo["enzymeIDs"]
			kcat = reactionInfo["kcat"][0]
			if reactionID in reactionEnzymesDict:
				for enzymeID in enzymeIDs:
					if enzymeID in reactionEnzymesDict[reactionID]:
						if kcat > reactionEnzymesDict[reactionID][enzymeID]:
							reactionEnzymesDict[reactionID][enzymeID] = kcat
		return reactionEnzymesDict

	def temperatureAdjustedKcat(self, reactionInfo):
		return 0


class ConcentrationUpdates(object):
	def __init__(self, concDict, equilibriumReactions, nutrientData):
		self.units = units.getUnit(concDict.values()[0])
		self.defaultConcentrationsDict = dict((key, concDict[key].asNumber(self.units)) for key in concDict)
		self.nutrientData = nutrientData

		self.moleculeScaleFactors = {
			"L-ALPHA-ALANINE[c]": 2.,
			"ARG[c]": 3.,
			"ASN[c]": 2.,
			"L-ASPARTATE[c]": 2.,
			"CYS[c]": 2.,
			"GLT[c]": 1.1,
			"GLN[c]": 2.,
			"GLY[c]": 2.,
			"HIS[c]": 2.,
			"ILE[c]": 2.,
			"LEU[c]": 2.,
			"LYS[c]": 2.,
			"MET[c]": 2.,
			"PHE[c]": 2.,
			"PRO[c]": 2.,
			"SER[c]": 2.,
			"THR[c]": 2.,
			"TRP[c]": 2.,
			"TYR[c]": 2.,
			"L-SELENOCYSTEINE[c]": 2.,
			"VAL[c]": 2.,
		}

		self.moleculeSetAmounts = self._addMoleculeAmountsBasedOnKd(equilibriumReactions)

	def concentrationsBasedOnNutrients(self, nutrientsLabel = None, nutrientsToInternalConc = None):
		concentrationsDict = self.defaultConcentrationsDict.copy()

		poolIds = sorted(concentrationsDict.keys())
		concentrations = self.units * np.array([concentrationsDict[k] for k in poolIds])

		if nutrientsLabel == None:
			return dict(zip(poolIds, concentrations))

		nutrientFluxes = {
			"importConstrainedExchangeMolecules": self.nutrientData["importConstrainedExchangeMolecules"][nutrientsLabel],
			"importUnconstrainedExchangeMolecules": self.nutrientData["importUnconstrainedExchangeMolecules"][nutrientsLabel],
		}

		concUpdates = None
		if nutrientsToInternalConc and nutrientsLabel in nutrientsToInternalConc:
			concUpdates = nutrientsToInternalConc[nutrientsLabel]
			concDict = dict(zip(poolIds, concentrations))
			concDict.update(concUpdates)
		else:
			for molecule, scaleFactor in self.moleculeScaleFactors.iteritems():
				if self._isNutrientExchangePresent(nutrientFluxes, molecule):
					concentrations[poolIds.index(molecule)] *= scaleFactor
			concDict = dict(zip(poolIds, concentrations))

		for moleculeName, setAmount in self.moleculeSetAmounts.iteritems():
			if concUpdates != None and moleculeName in concUpdates:
				continue
			if self._isNutrientExchangePresent(nutrientFluxes, moleculeName):
				concDict[moleculeName] = np.max((
					concDict.get(moleculeName, 0 * (units.mol / units.L)).asNumber(units.mol / units.L),
					setAmount.asNumber(units.mol / units.L)
					)) * (units.mol / units.L)

		return concDict

	def _isNutrientExchangePresent(self, nutrientFluxes, molecule):
		if molecule in nutrientFluxes["importUnconstrainedExchangeMolecules"]:
			return True

		if molecule in nutrientFluxes["importConstrainedExchangeMolecules"]:
			if nutrientFluxes["importConstrainedExchangeMolecules"][molecule].asNumber() > 0:
				return True

		return False

	def _addMoleculeAmountsBasedOnKd(self, equilibriumReactions):
		moleculeSetAmounts = {}
		for reaction in equilibriumReactions:
			# We only want to do this for species with standard Michaelis-Menten kinetics initially
			if len(reaction["stoichiometry"]) != 3:
				continue

			rev = reaction["reverse rate"]
			fwd = reaction["forward rate"]
			Kd = (units.mol / units.L) * (rev / fwd)

			moleculeName = [x["molecule"].encode("utf-8") for x in reaction["stoichiometry"] if x["type"] == "metabolite"][0]
			amountToSet = 0.
			if moleculeName in moleculeSetAmounts and moleculeSetAmounts[moleculeName] > Kd.asNumber(units.mol / units.L):
				amountToSet = moleculeSetAmounts[moleculeName]
			else:
				amountToSet = 1e-4#Kd.asNumber(units.mol / units.L)
			moleculeSetAmounts[moleculeName + "[p]"] = amountToSet * (units.mol / units.L)
			moleculeSetAmounts[moleculeName + "[c]"] = amountToSet * (units.mol / units.L)
		return moleculeSetAmounts
