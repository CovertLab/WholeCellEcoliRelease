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
from wholecell.utils.write_metabolic_constraints_file import writeMetabolicConstraintsFile
import wholecell
import os
import numpy as np
import collections
import warnings
import sympy as sp
from copy import copy

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

USE_ALL_CONSTRAINTS = False # False will remove problematic constraints from objective
CONSTRAINTS_TO_DISABLE = [
	"R601-RXN-FUM/REDUCED-MENAQUINONE//SUC/CPD-9728.38.",
	"SUCCINATE-DEHYDROGENASE-UBIQUINONE-RXN-SUC/UBIQUINONE-8//FUM/CPD-9956.31.",
	]

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

		# Create vector of metabolite target concentrations

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

		metaboliteIDs.append("PI[c]") # TODO: find a real number
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
		reactionCatalysts = {}
		catalystsList = []

		for reaction in raw_data.reactions:
			reactionID = reaction["reaction id"]
			stoich = reaction["stoichiometry"]
			reversible = reaction["is reversible"]

			if len(stoich) <= 1:
				raise Exception("Invalid biochemical reaction: {}, {}".format(reactionID, stoich))

			reactionStoich[reactionID] = stoich

			catalystsForThisRxn = []
			for catalyst in reaction["catalyzed by"]:
				try:
					catalystWithLoc = (catalyst + "[" + sim_data.getter.getLocation([catalyst])[0][0] + "]").encode("utf-8")
					catalystsForThisRxn.append(catalystWithLoc)
					catalystsList.append(catalystWithLoc)
				# If we don't have the catalyst in our reconstruction, drop it
				except KeyError:
					pass

			if len(catalystsForThisRxn) > 0:
				reactionCatalysts[reactionID] = catalystsForThisRxn

			# Add the reverse reaction
			if reversible:
				reverseReactionID = reverseReactionString.format(reactionID)
				reactionStoich[reverseReactionID] = {
					moleculeID:-stoichCoeff
					for moleculeID, stoichCoeff in reactionStoich[reactionID].viewitems()
					}

				reversibleReactions.append(reactionID)
				if len(catalystsForThisRxn) > 0:
					reactionCatalysts[reverseReactionID] = reactionCatalysts[reactionID]


		constraintDict = {}
		constraintIdList = []
		enzymeIdList = []
		kineticsSubstratesList = []
		reactionsToConstraintsDict = {}
		for constraint in raw_data.enzymeKinetics:
			if constraint["rateEquationType"] == "custom":
				continue

			constraintId = constraint["reactionID"].encode("utf-8") + "__" + constraint["enzymeIDs"].encode("utf-8") + "__%f" % (constraint["kcat"].asNumber(1 / units.s))
			assert len(constraint["Concentration Substrates"]) == len(constraint["kM"]) + len(constraint["kI"]), "Concentration Substrates are wrong length"
			assert constraintId not in constraintIdList, "constraintId already exists"

			# Get compartment for enzyme
			enzymeId = (constraint["enzymeIDs"] + "[" + sim_data.getter.getLocation([constraint["enzymeIDs"]])[0][0] + "]").encode("utf-8")
			constraint["enzymeIDs"] = enzymeId
			assert enzymeId in reactionCatalysts[constraint["reactionID"]], "%s is not a catalyst for %s according to FBA reconstruction" % (enzymeId, constraint["reactionID"])

			# Get compartments for Concentration Substrates
			concentrationSubstrates = []
			for substrate in constraint["Concentration Substrates"]:
				# In current implementation, anything with a concentration exists in the cytosol
				substrateWithCompartment = substrate.encode("utf-8") + "[c]"
				if substrateWithCompartment not in self.concDict:
					raise Exception, "Don't have concentration for %s" % substrateWithCompartment
				concentrationSubstrates.append(substrateWithCompartment)
				kineticsSubstratesList.append(substrateWithCompartment)
			constraint["Concentration Substrates"] = concentrationSubstrates

			# Get compartments for substrates
			substrates = []
			for substrate in constraint["substrateIDs"]:
				# In current implementation, anything with a concentration exists in the cytosol
				substrateFound = False
				for rxnSubstrate in reactionStoich[constraint["reactionID"]]:
					if rxnSubstrate.startswith(substrate + "["):
						substrateFound = True
						substrates.append(rxnSubstrate.encode("utf-8"))
				if not substrateFound:
					raise Exception, "Could not find compartment for substrate %s" % substrate
			constraint["substrateIDs"] = substrates

			# Adjust kcat for temperature
			temperature = constraint["Temp"]
			# If temperature not reported, assume 25 C
			if type(temperature) == str:
				temperature = 25
			constraint["kcatAdjusted"] = 2**((37. - temperature) / 10.) * constraint["kcat"]

			# Fix reactionID based on directionality
			stoichiometry = reactionStoich[constraint["reactionID"]]
			stoichVals = []
			for substrate in constraint["substrateIDs"]:
				stoichVals.append(stoichiometry[substrate])
			stoichVals = np.array(stoichVals)

			if np.all(stoichVals < 0):
				forward = True
			elif np.all(stoichVals > 0):
				forward = False
			else:
				raise Exception, "Have data for some reactants and some products (this is an inconsistency)"

			constraint["reactionID"] = constraint["reactionID"].encode("utf-8")
			if forward == False:
				constraint["reactionID"] = reverseReactionString.format(constraint["reactionID"])

			# Get rid of constraints for reverse reactions that the FBA reconstruction says should not exist
			# (i.e., if the FBA reconstruction says the reaction is irreversible but we have a constraint on 
			#  the reverse reaction, drop the constraint)
			if constraint["reactionID"] not in reactionStoich:
				continue

			enzymeIdList.append(enzymeId)
			constraintIdList.append(constraintId)
			constraintDict[constraintId] = constraint
			if constraint["reactionID"] not in reactionsToConstraintsDict:
				reactionsToConstraintsDict[constraint["reactionID"]] = []
			reactionsToConstraintsDict[constraint["reactionID"]].append(constraintId)

		constraintIdList = sorted(constraintIdList)
		constrainedReactionList = sorted(reactionsToConstraintsDict)
		kineticsSubstratesList = sorted(set(kineticsSubstratesList))
		enzymeIdList = sorted(set(enzymeIdList))

		# split out reactions that are kinetically constrained and that have more than one enzyme that catalyzes the reaction
		for rxn in constrainedReactionList:
			catalysts = reactionCatalysts[rxn]
			if len(catalysts) > 1:
				for catalyst in catalysts:
					if rxn.endswith(" (reverse)"):
						newReaction = reverseReactionString.format("%s__%s" % (rxn[:-10], catalyst[:-3]))
					else:
						newReaction = "%s__%s" % (rxn, catalyst[:-3])

					if rxn in reversibleReactions:
						reversibleReactions.append(newReaction)

					reactionStoich[newReaction] = copy(reactionStoich[rxn])
					reactionCatalysts[newReaction] = [catalyst]
					for constraint in reactionsToConstraintsDict[rxn]:
						if constraintDict[constraint]["enzymeIDs"] == catalyst:
							constraintDict[constraint]["reactionID"] = newReaction
							if newReaction not in reactionsToConstraintsDict:
								reactionsToConstraintsDict[newReaction] = []
							reactionsToConstraintsDict[newReaction].append(constraint)

				reactionStoich.pop(rxn)
				reactionCatalysts.pop(rxn)
				reactionsToConstraintsDict.pop(rxn)
				if rxn in reversibleReactions:
					reversibleReactions.pop(reversibleReactions.index(rxn))

		constrainedReactionList = sorted(reactionsToConstraintsDict)
		catalystsList = sorted(set(catalystsList))
		reactionCatalystsList = sorted(reactionCatalysts)

		# Create catalysis matrix (to be used in the simulation)
		catalysisMatrixI = []
		catalysisMatrixJ = []
		catalysisMatrixV = []

		for row, reaction in enumerate(reactionCatalystsList):
			for catalyst in reactionCatalysts[reaction]:
				col = catalystsList.index(catalyst)
				catalysisMatrixI.append(row)
				catalysisMatrixJ.append(col)
				catalysisMatrixV.append(1)

		catalysisMatrixI = np.array(catalysisMatrixI)
		catalysisMatrixJ = np.array(catalysisMatrixJ)
		catalysisMatrixV = np.array(catalysisMatrixV)

#		shape = (catalysisMatrixI.max() + 1, catalysisMatrixJ.max() + 1)
#		catalysisMatrix = np.zeros(shape, np.float64)
#		catalysisMatrix[catalysisMatrixI, catalysisMatrixJ] = catalysisMatrixV

		# Create constraint to reaction matrix (to be used in the simulation)
		constraintToReactionMatrixI = []
		constraintToReactionMatrixJ = []
		constraintToReactionMatrixV = []

		for row, reaction in enumerate(constrainedReactionList):
			for constraintId in reactionsToConstraintsDict[reaction]:
				col = constraintIdList.index(constraintId)
				constraintToReactionMatrixI.append(row)
				constraintToReactionMatrixJ.append(col)
				constraintToReactionMatrixV.append(1)

		constraintToReactionMatrixI = np.array(constraintToReactionMatrixI)
		constraintToReactionMatrixJ = np.array(constraintToReactionMatrixJ)
		constraintToReactionMatrixV = np.array(constraintToReactionMatrixV)

#		shape = (constraintToReactionMatrixI.max() + 1, constraintToReactionMatrixJ.max() + 1)
#		constraintToReactionMatrix = np.zeros(shape, np.float64)
#		constraintToReactionMatrix[constraintToReactionMatrixI, constraintToReactionMatrixJ] = constraintToReactionMatrixV

		# Use Sympy to create vector function that returns all kinetic constraints
		kineticsSubstrates = sp.symbols(["kineticsSubstrates[%d]" % idx for idx in xrange(len(kineticsSubstratesList))])
		enzymes = sp.symbols(["enzymes[%d]" % idx for idx in xrange(len(enzymeIdList))])
		constraints = [sp.symbol.S.Zero] * len(constraintIdList)
		constraintIsKcatOnly = np.zeros(len(constraintIdList))

		for constraintIdx, constraintId in enumerate(constraintIdList):
			constraint = constraintDict[constraintId]
			if constraint["rateEquationType"] == "custom":
				raise Exception

			enzymeIdx = enzymeIdList.index(constraint["enzymeIDs"])
			constraints[constraintIdx] = constraint["kcatAdjusted"].asNumber(1 / units.s) * enzymes[enzymeIdx]
			if len(constraint["kM"]) == 0 and len(constraint["kI"]) == 0:
				constraintIsKcatOnly[constraintIdx] = 1

			concSubstratePos = 0
			for kM in constraint["kM"]:
				kineticsSubstrateIdx = kineticsSubstratesList.index(constraint["Concentration Substrates"][concSubstratePos])
				S = kineticsSubstrates[kineticsSubstrateIdx]
				constraints[constraintIdx] *= (S / (S + kM))
				concSubstratePos += 1

			for kI in constraint["kI"]:
				kineticsSubstrateIdx = kineticsSubstratesList.index(constraint["Concentration Substrates"][concSubstratePos])
				I = kineticsSubstrates[kineticsSubstrateIdx]
				constraints[constraintIdx] *= (kI / (kI + I))
				concSubstratePos += 1


		constraints = sp.Matrix(constraints)

		constraintsFile = os.path.join(
			os.path.dirname(os.path.dirname(wholecell.__file__)),
			"reconstruction", "ecoli", "dataclasses", "process", "metabolism_constraints.py"
			)
		writeMetabolicConstraintsFile(constraintsFile, constraints)

		# Properties for FBA reconstruction
		self.reactionStoich = reactionStoich
		self.nutrientsTimeSeries = sim_data.nutrientsTimeSeries
		self.maintenanceReaction = {"ATP[c]": -1, "WATER[c]": -1, "ADP[c]": +1, "PI[c]": +1, "PROTON[c]": +1,}
		self.reversibleReactions = reversibleReactions

		# Properties for catalysis matrix (to set hard bounds)
		self.reactionCatalysts = reactionCatalysts
		self.catalystsList = catalystsList
		self.reactionCatalystsList = reactionCatalystsList
		self.catalysisMatrixI = catalysisMatrixI
		self.catalysisMatrixJ = catalysisMatrixJ
		self.catalysisMatrixV = catalysisMatrixV

		# Properties for setting flux targets
		self.constraintIdList = constraintIdList
		self.constrainedReactionList = constrainedReactionList
		self.constraintToReactionMatrixI = constraintToReactionMatrixI
		self.constraintToReactionMatrixJ = constraintToReactionMatrixJ
		self.constraintToReactionMatrixV = constraintToReactionMatrixV
		self.enzymeIdList = enzymeIdList
		self.kineticsSubstratesList = kineticsSubstratesList
		self.constraintDict = constraintDict
		self.reactionsToConstraintsDict = reactionsToConstraintsDict
		self.constraintIsKcatOnly = constraintIsKcatOnly
		self.useAllConstraints = USE_ALL_CONSTRAINTS
		self.constraintsToDisable = CONSTRAINTS_TO_DISABLE



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


class ConcentrationUpdates(object):
	def __init__(self, concDict, equilibriumReactions, nutrientData):
		self.units = units.getUnit(concDict.values()[0])
		self.defaultConcentrationsDict = dict((key, concDict[key].asNumber(self.units)) for key in concDict)
		self.nutrientData = nutrientData

		self.moleculeScaleFactors = {
			"L-ALPHA-ALANINE[c]": 2.,
			"ARG[c]": 2.,
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

		self.moleculeSetAmounts = self._addMoleculeAmounts(equilibriumReactions, self.defaultConcentrationsDict)

	def concentrationsBasedOnNutrients(self, nutrientsLabel = None, nutrientsToInternalConc = None):
		concentrationsDict = self.defaultConcentrationsDict.copy()

		metaboliteTargetIds = sorted(concentrationsDict.keys())
		concentrations = self.units * np.array([concentrationsDict[k] for k in metaboliteTargetIds])

		if nutrientsLabel == None:
			return dict(zip(metaboliteTargetIds, concentrations))

		nutrientFluxes = {
			"importConstrainedExchangeMolecules": self.nutrientData["importConstrainedExchangeMolecules"][nutrientsLabel],
			"importUnconstrainedExchangeMolecules": self.nutrientData["importUnconstrainedExchangeMolecules"][nutrientsLabel],
		}

		concDict = dict(zip(metaboliteTargetIds, concentrations))

		for moleculeName, setAmount in self.moleculeSetAmounts.iteritems():
			if self._isNutrientExchangePresent(nutrientFluxes, moleculeName) and (moleculeName[:-3] + "[c]" not in self.moleculeScaleFactors or moleculeName == "L-SELENOCYSTEINE[c]"):
				concDict[moleculeName] = setAmount
			if moleculeName in self.moleculeScaleFactors and self._isNutrientExchangePresent(nutrientFluxes, moleculeName[:-3] + "[p]"):
				concDict[moleculeName] = setAmount

		return concDict

	def _isNutrientExchangePresent(self, nutrientFluxes, molecule):
		if molecule in nutrientFluxes["importUnconstrainedExchangeMolecules"]:
			return True

		if molecule in nutrientFluxes["importConstrainedExchangeMolecules"]:
			if nutrientFluxes["importConstrainedExchangeMolecules"][molecule].asNumber() > 0:
				return True

		return False

	def _addMoleculeAmounts(self, equilibriumReactions, concDict):
		moleculeSetAmounts = {}
		for reaction in equilibriumReactions:
			# We only want to do this for species with standard Michaelis-Menten kinetics initially
			if len(reaction["stoichiometry"]) != 3:
				continue

			moleculeName = [x["molecule"].encode("utf-8") for x in reaction["stoichiometry"] if x["type"] == "metabolite"][0]
			amountToSet = 1e-4
			moleculeSetAmounts[moleculeName + "[p]"] = amountToSet * self.units
			moleculeSetAmounts[moleculeName + "[c]"] = amountToSet * self.units

		for moleculeName, scaleFactor in self.moleculeScaleFactors.iteritems():
			moleculeSetAmounts[moleculeName] = scaleFactor * concDict[moleculeName] * self.units
		return moleculeSetAmounts
