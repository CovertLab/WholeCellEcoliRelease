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

ILE_LEU_CONCENTRATION = 3.0e-4 # mmol/L
ILE_FRACTION = 0.360 # the fraction of iso/leucine that is isoleucine; computed from our monomer data
ECOLI_PH = 7.2

PPI_CONCENTRATION = 0.5e-3 # M, multiple sources

EXCHANGE_UNITS = units.mmol / units.g / units.h

class Metabolism(object):
	""" Metabolism """

	def __init__(self, raw_data, sim_data):
		self._buildBiomass(raw_data, sim_data)
		self._buildMetabolism(raw_data, sim_data)


	def _buildBiomass(self, raw_data, sim_data):
		wildtypeIDs = set(entry["molecule id"] for entry in raw_data.biomass)
		# TODO: unjank this

		# Load the biomass function flat file as a dict
		biomassFunction = {entry['molecule id']:entry['coefficient'] for entry in raw_data.biomass}

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

		self.biomassFunction = biomassFunction
		self.concentrationUpdates = ConcentrationUpdates(dict(zip(
			metaboliteIDs,
			(units.mol / units.L) * np.array(metaboliteConcentrations)
			)))
		envFirstTimePoint = sim_data.envDict[sim_data.environment][0][-1]
		self.concDict = self.concentrationUpdates.concentrationsBasedOnNutrients(envFirstTimePoint)

	def _buildMetabolism(self, raw_data, sim_data):
		# Build the matrices/vectors for metabolism (FBA)

		# These may be modified/extended later, but should provide the basic
		# data structures

		reactionStoich = {}
		reversibleReactions = []
		reactionEnzymes = {}
		reactionRates = {}

		validEnzymeIDs = set([])
		validProteinIDs = ['{}[{}]'.format(x['id'],location) for x in raw_data.proteins for location in x['location']]
		validProteinComplexIDs = ['{}[{}]'.format(x['id'],location) for x in raw_data.proteinComplexes for location in x['location']]
		validEnzymeIDs.update(validProteinIDs)
		validEnzymeIDs.update(validProteinComplexIDs)
		validEnzymeCompartments = collections.defaultdict(set)

		for enzymeID in validEnzymeIDs:
			enzyme = enzymeID[:enzymeID.index("[")]
			location = enzymeID[enzymeID.index("[")+1:enzymeID.index("[")+2]

			validEnzymeCompartments[enzyme].add(location)

		for reaction in raw_data.reactions:
			reactionID = reaction["reaction id"]
			stoich = reaction["stoichiometry"]

			if len(stoich) <= 1:
				raise Exception("Invalid biochemical reaction: {}, {}".format(reactionID, stoich))

			reactionStoich[reactionID] = stoich

			# Assign reversibilty
			if reaction["is reversible"]:
				reversibleReactions.append(reactionID)

		reactionRateInfo = {}
		constraintIDs = []
		constraintToReactionDict = {}
		activeConstraintsDict = {}
		enzymesWithKineticInfo = set()
		enzymesWithKineticInfoDict = {}

		# Enzyme kinetics data
		for reaction in raw_data.enzymeKinetics:
			constraintID = reaction["constraintID"]
			constraintIDs.append(constraintID)

			constraintToReactionDict[constraintID] = reaction["reactionID"]
			activeConstraintsDict[constraintID] = reaction["constraintActive"]

			# If the enzymes don't already have a compartment tag, add one from the valid compartment list or [c] (cytosol) as a default
			new_reaction_enzymes = []
			for reactionEnzyme in reaction["enzymeIDs"]:
				if reactionEnzyme[-3:-2] !='[':
					if len(validEnzymeCompartments[reactionEnzyme]) > 0:
						new_reaction_enzymes.append(reactionEnzyme +'['+str(validEnzymeCompartments[reactionEnzyme].pop())+']')
					else:
						new_reaction_enzymes.append(reactionEnzyme + '[c]')
				else:
					new_reaction_enzymes.append(reactionEnzyme)

			reaction["enzymeIDs"] = new_reaction_enzymes

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


			reactionRateInfo[constraintID] = reaction
			for enzymeID in reaction["enzymeIDs"]:
				enzymesWithKineticInfo.add(enzymeID)


		enzymesWithKineticInfoDict["enzymes"] = list(enzymesWithKineticInfo)

		self.reactionStoich = reactionStoich
		self.envDict = sim_data.envDict
		self.reversibleReactions = reversibleReactions
		self.reactionRateInfo = reactionRateInfo
		self.enzymesWithKineticInfo = enzymesWithKineticInfoDict
		self.constraintIDs = constraintIDs
		self.constraintToReactionDict = constraintToReactionDict
		self.activeConstraintsDict = activeConstraintsDict

	def exchangeConstraints(self, exchangeIDs, coefficient, targetUnits, environment, time):
		newObjective = None
		while len(self.envDict[environment]) and time > self.envDict[environment][0][0]:
			self._unconstrainedExchangeMolecules = self.envDict[environment][0][1]["unconstrainedExchangeMolecules"]
			self._constrainedExchangeMolecules = self.envDict[environment][0][1]["constrainedExchangeMolecules"]
			concDict = self.concentrationUpdates.concentrationsBasedOnNutrients(self.envDict[environment][0][-1])
			newObjective = dict((key, concDict[key].asNumber(targetUnits)) for key in concDict)
			self.envDict[environment].popleft()

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

class ConcentrationUpdates(object):
	def __init__(self, concDict):
		self.units = units.getUnit(concDict.values()[0])
		self.defaultConcentrationsDict = dict((key, concDict[key].asNumber(self.units)) for key in concDict)

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

	def concentrationsBasedOnNutrients(self, nutrientFluxes = None):
		concentrationsDict = self.defaultConcentrationsDict.copy()

		poolIds = sorted(concentrationsDict.keys())
		concentrations = self.units * np.array([concentrationsDict[k] for k in poolIds])

		if nutrientFluxes == None:
			return dict(zip(poolIds, concentrations))

		for molecule, scaleFactor in self.moleculeScaleFactors.iteritems():
			if self._isNutrientExchangePresent(nutrientFluxes, molecule):
				concentrations[poolIds.index(molecule)] *= scaleFactor

		return dict(zip(poolIds, concentrations))

	def _isNutrientExchangePresent(self, nutrientFluxes, molecule):
		if molecule in nutrientFluxes["unconstrainedExchangeMolecules"]:
			return True

		if molecule in nutrientFluxes["constrainedExchangeMolecules"]:
			if nutrientFluxes["constrainedExchangeMolecules"][molecule].asNumber() > 0:
				return True

		return False