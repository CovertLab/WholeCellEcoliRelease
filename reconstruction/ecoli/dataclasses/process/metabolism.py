"""
SimulationData for metabolism process

TODO:
- improved estimate of ILE/LEU abundance or some external data point
- implement L1-norm minimization for AA concentrations
- find concentration for PI[c]
- add (d)NTP byproduct concentrations
- include custom kinetic constraints when reading from raw_data.enzymeKinetics

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 03/06/2015
"""

from __future__ import division

from wholecell.utils import units
import wholecell
import os
import numpy as np
import sympy as sp
from copy import copy

ILE_LEU_CONCENTRATION = 3.0e-4 # mmol/L
ILE_FRACTION = 0.360 # the fraction of iso/leucine that is isoleucine; computed from our monomer data
ECOLI_PH = 7.2

PPI_CONCENTRATION = 0.5e-3 # M, multiple sources

USE_ALL_CONSTRAINTS = False # False will remove defined constraints from objective

reverseReactionString = "{} (reverse)"

# threshold (units.mmol / units.L) separates concentrations that are import constrained with
# max flux = 0 from unconstrained molecules.
# TODO (Eran) remove this once a transport kinetics process is operating
IMPORT_CONSTRAINT_THRESHOLD =  0.01

class Metabolism(object):
	""" Metabolism """

	def __init__(self, raw_data, sim_data):
		# set solver and kinetic objective weight
		self.solver = "glpk-linear"
		if "linear" in self.solver:
			self.kinetic_objective_weight = sim_data.constants.metabolismKineticObjectiveWeightLinear
		else:
			self.kinetic_objective_weight = sim_data.constants.metabolismKineticObjectiveWeightQuadratic

		# lists of molecules whose presence modifies glc's upper bound for FBA import constraint, whose default is 20 (mmol/g DCW/hr).
		# This is implemented to reproduce glc maximum uptake previous instantiated in environment files, but now done explicitly here.
		# TODO (Eran) transport should do away with these conditional requirements for determining GLC flux.
		self.glc_vmax_conditions = [
			# if any of these molecules are ABSENT, GLC upper bound is set to 0 (mmol/g DCW/hr)
			['GLC[p]'],
			# if any of these molecules are ABSENT, GLC upper bound is set to 100 (mmol/g DCW/hr)
			['OXYGEN-MOLECULE[p]'],
			# if any of these molecules are ABSENT, GLC upper bound is set to 10 (mmol/g DCW/hr)
			['CA+2[p]', 'MG+2[p]', 'PI[p]'],
			# if any of these molecules are PRESENT, GLC upper bound is set to 10 (mmol/g DCW/hr)
			['CPD-183[p]', 'INDOLE[p]', 'NITRATE[p]', 'NITRITE[p]', 'CPD-520[p]', 'TUNGSTATE[p]'],
		]

		self.all_external_exchange_molecules = self._getAllExternalExchangeMolecules(raw_data)
		self._secretion_exchange_molecules = self._getSecretionExchangeMolecules(raw_data)
		self._exchange_data_dict = self._getExchangeDataDict(raw_data, sim_data)
		self.import_exchange =[]
		self.import_constraint = []

		# initialize exchange_data, import_exchange, and import_constraints
		self.exchange_data = self.getExchangeData('minimal')
		self.saveImportConstraints(self.exchange_data)

		self._buildBiomass(raw_data, sim_data)
		self._buildMetabolism(raw_data, sim_data)


	def saveImportConstraints(self, exchange_data):
		'''
		Saves import_constraint and import_exchange for the fba_results listener.
		import_constraint saves all importConstrainedExchangeMolecules as true, and the rest as false.
		import_exchange saves all importExchangeMolecules as true, and the rest as false.
		'''

		# molecules from all_external_exchange_molecules set to 'true' if they are current importExchangeMolecules.
		self.import_exchange = [
			molecule_id in exchange_data['importExchangeMolecules']
			for molecule_id in self.all_external_exchange_molecules
			]

		# molecules from all_external_exchange_molecules set to 'true' if they are current importConstrainedExchangeMolecules.
		self.import_constraint = [
			molecule_id in exchange_data['importConstrainedExchangeMolecules']
			for molecule_id in self.all_external_exchange_molecules
			]


	def exchangeDataFromConcentrations(self, molecules):
		'''
		Update importExchangeMolecules for FBA based on current nutrient concentrations.
		This provides a simple type of transport to accommodate changing nutrient
		concentrations in the environment. Transport is modeled as a binary switch:
		When there is a high concentrations of environment nutrients, transporters
		are unconstrained and nutrients are transported as needed by metabolism.
		When concentrations fall below the threshold, that nutrient's transport
		is constrained to max flux of 0.

		Five categories of molecules set up the FBA problem space. These include:
			- externalExchangeMolecules: All exchange molecules, both import and secretion exchanged molecules.
			- importExchangeMolecules: molecules that can be imported from the environment into the cell.
			- importConstrainedExchangeMolecules: exchange molecules that have an upper bound on their flux.
			- importUnconstrainedExchangeMolecules: exchange molecules that do not have an upper bound on their flux.
			- secretionExchangeMolecules: molecules that can be secreted by the	cell into the environment.
		'''

		externalExchangeMolecules = set()
		importExchangeMolecules = set()
		importUnconstrainedExchangeMolecules = []
		importConstrainedExchangeMolecules = {}
		secretionExchangeMolecules = self._secretion_exchange_molecules

		#remove molecules with low concentration
		nonzero_molecules = {molecule_id:concentration for molecule_id, concentration in molecules.items() if concentration >= 0.00001}

		for molecule_id, concentration in nonzero_molecules.iteritems():

			# skip if concentration is 0. Do not include these nonexistent molecules in FBA problem definition.
			if molecule_id != 'GLC[p]' and concentration == 0:
				continue

			elif concentration < IMPORT_CONSTRAINT_THRESHOLD:
				importConstrainedExchangeMolecules[molecule_id] = 0 * (units.mmol / units.g / units.h)

			# The logic below is used to change GLC's upper bound flux based on what nutrients are present in the environment.
			# The order of this logic is important. First, if GLC is absent from the environment, the upper bound flux is 0.
			# If molecules in condition[1] (oxygen) are absent, the upper bound flux is 100 to match anaerobic glucose uptake.
			# If either molecules in condition[2] are absent or condition[3] are present, the upper bound flux is 10.
			# Finally, if none of these conditions are true, the upper bound is set to the default of 20.
			elif molecule_id == 'GLC[p]':
				# if any molecule in glc_vmax_conditions[0] is ABSENT:
				if not all(molecule in nonzero_molecules for molecule in self.glc_vmax_conditions[0]):
					importConstrainedExchangeMolecules[molecule_id] = 0 * (units.mmol / units.g / units.h)
				# if any molecule in glc_vmax_conditions[1] is ABSENT:
				elif not all(molecule in nonzero_molecules for molecule in self.glc_vmax_conditions[1]):
					importConstrainedExchangeMolecules[molecule_id] = 100 * (units.mmol / units.g / units.h)
				# if any molecule in glc_vmax_conditions[2] is ABSENT:
				elif not all(molecule in nonzero_molecules for molecule in self.glc_vmax_conditions[2]):
					importConstrainedExchangeMolecules[molecule_id] = 10 * (units.mmol / units.g / units.h)
				# if any molecule in glc_vmax_conditions[3] is PRESENT:
				elif any(molecule in nonzero_molecules for molecule in self.glc_vmax_conditions[3]):
					importConstrainedExchangeMolecules[molecule_id] = 10 * (units.mmol / units.g / units.h)
				else:
					importConstrainedExchangeMolecules[molecule_id] = 20 * (units.mmol / units.g / units.h)

			# add to import unconstrained if concentration >= threshold and not GLC
			else:
				importUnconstrainedExchangeMolecules.append(molecule_id)

			importExchangeMolecules.add(molecule_id)
			externalExchangeMolecules.add(molecule_id)

		for molecule_id in secretionExchangeMolecules:
			externalExchangeMolecules.add(molecule_id)

		return {
			"externalExchangeMolecules": list(externalExchangeMolecules),
			"importExchangeMolecules": list(importExchangeMolecules),
			"importConstrainedExchangeMolecules": importConstrainedExchangeMolecules,
			"importUnconstrainedExchangeMolecules": importUnconstrainedExchangeMolecules,
			"secretionExchangeMolecules": secretionExchangeMolecules,
		}


	def _getExchangeDataDict(self, raw_data, sim_data):
		'''
		Returns a dictionary of exchange_data for the initial environment listed in condition.media. This dictionary
		is used to quickly pull up  exchange data for these different environments by their name.
		'''

		self.environment_dict = sim_data.external_state.environment.environment_dict

		externalExchangeMolecules = {}
		importExchangeMolecules = {}
		importConstrainedExchangeMolecules = {}
		importUnconstrainedExchangeMolecules = {}
		secretionExchangeMolecules = self._secretion_exchange_molecules

		for environment_name, molecules in self.environment_dict.iteritems():

			exchange_data = self.exchangeDataFromConcentrations(molecules)

			externalExchangeMolecules[environment_name] = exchange_data['externalExchangeMolecules']
			importExchangeMolecules[environment_name] = exchange_data['importExchangeMolecules']
			importConstrainedExchangeMolecules[environment_name] = exchange_data['importConstrainedExchangeMolecules']
			importUnconstrainedExchangeMolecules[environment_name] = exchange_data['importUnconstrainedExchangeMolecules']

		return {
			"externalExchangeMolecules": externalExchangeMolecules,
			"importExchangeMolecules": importExchangeMolecules,
			"importConstrainedExchangeMolecules": importConstrainedExchangeMolecules,
			"importUnconstrainedExchangeMolecules": importUnconstrainedExchangeMolecules,
			"secretionExchangeMolecules": secretionExchangeMolecules,
		}


	def getExchangeData(self, environment_label):
		'''
		Returns exchange_data for a given environment_label saved in _exchange_data_dict.
		'''

		externalExchangeMolecules = self._exchange_data_dict['externalExchangeMolecules'][environment_label]
		importExchangeMolecules = self._exchange_data_dict['importExchangeMolecules'][environment_label]
		importConstrainedExchangeMolecules = self._exchange_data_dict['importConstrainedExchangeMolecules'][environment_label]
		importUnconstrainedExchangeMolecules = self._exchange_data_dict['importUnconstrainedExchangeMolecules'][environment_label]
		secretionExchangeMolecules = self._exchange_data_dict['secretionExchangeMolecules']

		return {
			"externalExchangeMolecules": externalExchangeMolecules,
			"importExchangeMolecules": importExchangeMolecules,
			"importConstrainedExchangeMolecules": importConstrainedExchangeMolecules,
			"importUnconstrainedExchangeMolecules": importUnconstrainedExchangeMolecules,
			"secretionExchangeMolecules": secretionExchangeMolecules,
		}


	def _getSecretionExchangeMolecules(self, raw_data):
		'''
		get list of all secretion exchange molecules from raw data
		'''
		secretionExchangeMolecules = []
		for secretion in raw_data.secretions:
			if secretion["lower bound"] and secretion["upper bound"]:
				# "non-growth associated maintenance", not included in our metabolic model
				continue
			else:
				secretionExchangeMolecules.append(secretion["molecule id"])

		return secretionExchangeMolecules


	def _getAllExternalExchangeMolecules(self, raw_data):
		'''
		get list of all external exchange molecules from raw data
		'''
		externalExchangeData = []
		# initiate all molecules with 0 concentrations
		for row in raw_data.condition.environment_molecules:
			externalExchangeData.append(row["molecule id"])

		return externalExchangeData


	def _buildBiomass(self, raw_data, sim_data):
		wildtypeIDs = set(entry["molecule id"] for entry in raw_data.biomass)

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
		ileRelative = ILE_FRACTION
		leuRelative = 1 - ileRelative

		metaboliteIDs.append("ILE[c]")
		metaboliteConcentrations.append(ileRelative * ILE_LEU_CONCENTRATION)

		metaboliteIDs.append("LEU[c]")
		metaboliteConcentrations.append(leuRelative * ILE_LEU_CONCENTRATION)

		# CYS/SEC/GLY: concentration based on other amino acids
		aaConcentrations = []

		for aaIndex, aaID in enumerate(sim_data.amino_acid_1_to_3_ordered.values()):
			if aaID in metaboliteIDs:
				metIndex = metaboliteIDs.index(aaID)
				aaConcentrations.append(metaboliteConcentrations[metIndex])

		aaSmallestConc = min(aaConcentrations)

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

		for dntpIndex, dntpID in enumerate(sim_data.moleculeGroups.dNtpIds):
			if dntpID in metaboliteIDs:
				metIndex = metaboliteIDs.index(dntpID)
				dntpConcentrations.append(metaboliteConcentrations[metIndex])

		dntpSmallestConc = min(dntpConcentrations)

		metaboliteIDs.append("DGTP[c]")
		metaboliteConcentrations.append(dntpSmallestConc)

		# H: from reported pH
		hydrogenConcentration = 10**(-ECOLI_PH)

		metaboliteIDs.append("PROTON[c]")
		metaboliteConcentrations.append(hydrogenConcentration)

		# PPI: multiple sources report 0.5 mM
		metaboliteIDs.append("PPI[c]")
		metaboliteConcentrations.append(PPI_CONCENTRATION)

		metaboliteIDs.append("PI[c]")
		metaboliteConcentrations.append(PPI_CONCENTRATION)

		# Add byproducts with no annotated concentration to force recycling
		metaboliteIDs.append("UMP[c]")
		metaboliteConcentrations.append(2.40e-5)

		# include metabolites that are part of biomass
		for key, value in sim_data.mass.getBiomassAsConcentrations(sim_data.doubling_time).iteritems():
			metaboliteIDs.append(key)
			metaboliteConcentrations.append(value.asNumber(units.mol / units.L))

		# save concentrations as class variables
		self.concentrationUpdates = ConcentrationUpdates(dict(zip(
			metaboliteIDs,
			(units.mol / units.L) * np.array(metaboliteConcentrations)
			)),
			raw_data.equilibriumReactions,
			self._exchange_data_dict,
		)
		self.concDict = self.concentrationUpdates.concentrationsBasedOnNutrients("minimal")
		self.nutrientsToInternalConc = {}
		self.nutrientsToInternalConc["minimal"] = self.concDict.copy()

	def _buildMetabolism(self, raw_data, sim_data):
		"""
		Build the matrices/vectors for metabolism (FBA)
		Reads in and stores reaction and kinetic constraint information
		"""

		# Initialize variables to store reaction information
		reactionStoich = {}			# dict with reactions as keys and dict with reaction stoich as values 
		reversibleReactions = []
		reactionCatalysts = {}		# dict with reactions as keys and list of catalysts as values
		catalystsList = []

		# Load and parse reaction information from raw_data
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

		# Initialize variables to store kinetic constraint information
		constraintDict = {}
		constraintIdList = []
		enzymeIdList = []
		kineticsSubstratesList = []
		reactionsToConstraintsDict = {}

		# Load and parse kinetic constraint information from raw_data
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
					# create new reaction name with enzyme appended to the end
					if rxn.endswith(" (reverse)"):
						newReaction = reverseReactionString.format("%s__%s" % (rxn[:-10], catalyst[:-3]))
					else:
						newReaction = "%s__%s" % (rxn, catalyst[:-3])

					# add the new reaction to appropriate lists and dicts
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

				# remove old reaction name
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

		self._kineticConstraints = str(constraints)
		self._compiledConstraints = None

		# Properties for FBA reconstruction
		self.reactionStoich = reactionStoich
		self.nutrientsTimeSeries = sim_data.external_state.environment.nutrients_time_series
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
		self.constraintsToDisable = [rxn["disabled reaction"] for rxn in raw_data.disabledKineticReactions]

	def getKineticConstraints(self, enzymes, substrates):
		'''
		Allows for dynamic code generation for kinetic constraint calculation
		for use in Metabolism process. Inputs should be unitless but the order
		of magnitude should match the kinetics parameters (umol/L/s).

		If trying to pickle sim_data object after function has been called,
		_compiledConstraints might not be able to be pickled.  See
		__getstate__(), __setstate__() comments on PR 111 to address.

		Returns np.array of floats of the kinetic constraint target for each
		reaction with kinetic parameters
		Inputs:
			enzymes (np.array of floats) - concentrations of enzymes associated
				with kinetics constraints
			substrates (np.array of floats) - concentrations of substrates
				associated with kinetics constraints
		'''

		if self._compiledConstraints is None:
			self._compiledConstraints = eval('lambda enzymes, kineticsSubstrates: np.array(%s)\n'
				% self._kineticConstraints, {'np': np}, {}
				)

		return self._compiledConstraints(enzymes, substrates)

	def exchangeConstraints(self, exchangeIDs, coefficient, targetUnits, currentNutrients, exchange_data, concModificationsBasedOnCondition = None):
		"""
		Called during Metabolism process
		Returns the homeostatic objective concentrations based on the current nutrients
		Returns levels for external molecules available to exchange based on the current nutrients
		"""

		newObjective = None

		self._unconstrainedExchangeMolecules = exchange_data["importUnconstrainedExchangeMolecules"]
		self._constrainedExchangeMolecules = exchange_data["importConstrainedExchangeMolecules"]

		concDict = self.concentrationUpdates.concentrationsBasedOnNutrients(currentNutrients, self.nutrientsToInternalConc)
		if concModificationsBasedOnCondition is not None:
			concDict.update(concModificationsBasedOnCondition)

		conversion = targetUnits.asUnit(self.concentrationUpdates.units)
		newObjective = dict((key, (val / conversion).asNumber()) for key, val in concDict.iteritems())

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

# Class used to update metabolite concentrations based on the current nutrient conditions
class ConcentrationUpdates(object):
	def __init__(self, concDict, equilibriumReactions, exchange_data_dict):
		self.units = units.getUnit(concDict.values()[0])
		self.defaultConcentrationsDict = dict((key, concDict[key].asNumber(self.units)) for key in concDict)
		self._exchange_data_dict = exchange_data_dict

		# factor of internal amino acid increase if maino acids present in nutrients
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

	# return adjustments to concDict based on nutrient conditions
	def concentrationsBasedOnNutrients(self, nutrientsLabel = None, nutrientsToInternalConc = None):
		concentrationsDict = self.defaultConcentrationsDict.copy()

		metaboliteTargetIds = sorted(concentrationsDict.keys())
		concentrations = self.units * np.array([concentrationsDict[k] for k in metaboliteTargetIds])

		if nutrientsLabel == None:
			return dict(zip(metaboliteTargetIds, concentrations))

		nutrientFluxes = {
			"importConstrainedExchangeMolecules": self._exchange_data_dict["importConstrainedExchangeMolecules"][nutrientsLabel],
			"importUnconstrainedExchangeMolecules": self._exchange_data_dict["importUnconstrainedExchangeMolecules"][nutrientsLabel],
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
