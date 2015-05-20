"""
SimulationData for metabolism process

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 03/06/2015
"""

from __future__ import division

from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray
import numpy as np
import collections

CELL_DENSITY = 1.1e3 # g/L
# from Baldwin WW, Myer R, Powell N, Anderson E, Koch AL. Buoyant
# density of Escherichia coli is determined solely by the osmolarity
# of the culture medium. Arch Microbiol. 1995 Aug164(2):155-7 p.156
# fig.3 & fig.2 (retrieved from Bionumbers)

METABOLITE_CONCENTRATIONS = { # mol / L # TODO: move to SQL # TODO: Add units!
	"GLT": 9.60e-2,
	"GLUTATHIONE": 1.70e-2,
	"FRUCTOSE-16-DIPHOSPHATE": 1.50e-2,
	"ATP": 9.60e-3,
	"UDP-N-ACETYL-D-GLUCOSAMINE": 9.20e-3,
	"UTP": 8.30e-3,
	"GTP": 4.90e-3,
	"TTP": 4.60e-3,
	"L-ASPARTATE": 4.20e-3,
	"VAL": 4.00e-3,
	"CPD-2961": 3.80e-3, # D-gluconate 6-phosphate
	"GLN": 3.80e-3,
	"CTP": 2.70e-3,
	"ALA": 2.60e-3,
	"NAD": 2.60e-3,
	"CPD-12575": 2.50e-3, # UDP-glucose
	"URIDINE": 2.10e-3,
	"CIT": 2.00e-3,
	"UDP": 1.80e-3,
	"MAL": 1.70e-3, # this is (S)-malate, but there is also (R)-malate
	"G3P": 1.50e-3,
	"L-CITRULLINE": 1.40e-3,
	"CO-A": 1.40e-3,
	"GLYCERATE": 1.40e-3,
	"CPD-13469": 1.20e-3, # alpha-D-glucosamine-6-phosphate
	"ACETYL-P": 1.10e-3,
	"GLC-D-LACTONE": 1.00e-3,
	"GDP": 6.80e-4,
	"ACETYL-COA": 6.10e-4,
	"CARBAMYUL-L-ASPARTATE": 5.90e-4,
	"ARG": 5.70e-4,
	"SUC": 5.70e-4,
	"UDP-GLUCURONATE": 5.70e-4,
	"ADP": 5.60e-4,
	"ASN": 5.10e-4,
	"2-KETOGLUTARATE": 4.40e-4,
	"LYS": 4.10e-4,
	"PRO": 3.90e-4,
	"DTDP": 3.80e-4,
	"DIHYDROXY-ACETONE-PHOSPHATE": 3.70e-4,
	"HOMO-CYS": 3.70e-4,
	"CMP": 3.60e-4,
	"AMP": 2.80e-4,
	"SUC-COA": 2.30e-4,
	"GUANINE": 1.90e-4,
	"PHOSPHO-ENOL-PYRUVATE": 1.80e-4,
	"S-ADENOSYLMETHIONINE": 1.80e-4,
	"THR": 1.80e-4,
	"FAD": 1.70e-4,
	"MET": 1.50e-4,
	"2-3-DIHYDROXYBENZOATE": 1.40e-4,
	"FUM": 1.20e-4,
	"NADPH": 1.20e-4,
	"PHENYL-PYRUVATE": 9.00e-5,
	"NADH": 8.30e-5,
	"N-ACETYL-D-GLUCOSAMINE-1-P": 8.20e-5,
	"HIS": 6.80e-5,
	"SER": 6.80e-5,
	"4-hydroxybenzoate": 5.20e-5,
	"DGMP": 5.10e-5,
	"GLYCEROL-3P": 4.90e-5,
	"N-ALPHA-ACETYLORNITHINE": 4.30e-5,
	"GLUCONATE": 4.20e-5,
	# "CAMP": 3.50e-5,
	"DCTP": 3.50e-5,
	"MALONYL-COA": 3.50e-5,
	"TYR": 2.90e-5,
	"GMP": 2.40e-5,
	"ACETOACETYL-COA": 2.20e-5,
	"RIBOFLAVIN": 1.90e-5,
	"PHE": 1.80e-5,
	"CIS-ACONITATE": 1.60e-5,
	"DATP": 1.60e-5,
	"CYTOSINE": 1.40e-5,
	"SHIKIMATE": 1.40e-5,
	"HISTIDINOL": 1.30e-5,
	"DI-H-OROTATE": 1.20e-5,
	"QUINOLINATE": 1.20e-5,
	"TRP": 1.20e-5,
	"L-ORNITHINE": 1.00e-5,
	"DAMP": 8.80e-6,
	"APS": 6.60e-6,
	# "MYO-INOSITOL": 5.70e-6, # can't be formed by the reaction network
	"PROPIONYL-COA": 5.30e-6,
	"ADP-D-GLUCOSE": 4.30e-6,
	"ANTHRANILATE": 3.50e-6,
	"DEOXYADENOSINE": 2.80e-6,
	"CYTIDINE": 2.60e-6,
	"NADP": 2.10e-6,
	"GUANOSINE": 1.60e-6,
	"ADENINE": 1.50e-6,
	"DEOXYGUANOSINE": 5.20e-7,
	"ADENOSINE": 1.30e-7,
	}

ILE_LEU_CONCENTRATION = 3.0e-4 # mmol/L
ILE_FRACTION = 0.360 # the fraction of iso/leucine that is isoleucine; computed from our monomer data
ECOLI_PH = 7.2

PPI_CONCENTRATION = 0.5e-3 # M, multiple sources

class Metabolism(object):
	""" Metabolism """

	def __init__(self, raw_data, sim_data):
		self._buildBiomass(raw_data, sim_data)
		self._buildMetabolism(raw_data, sim_data)


	def _buildBiomass(self, raw_data, sim_data):
		wildtypeIDs = set(entry["molecule id"] for entry in raw_data.biomass)
		# TODO: unjank this

		# Create vector of metabolite pools (concentrations)

		# Since the data only covers certain metabolites, we need to rationally
		# expand the dataset to include the other molecules in the biomass
		# function.

		# First, load in metabolites that do have concentrations, then assign
		# compartments according to those given in the biomass objective.  Or,
		# if there is no compartment, assign it to the cytoplasm.

		metaboliteIDs = []
		metaboliteConcentrations = []

		wildtypeIDtoCompartment = {
			wildtypeID[:-3] : wildtypeID[-3:]
			for wildtypeID in wildtypeIDs
			} # this assumes biomass reaction components only exist in a single compartment

		for metaboliteID, concentration in METABOLITE_CONCENTRATIONS.viewitems():
			if metaboliteID in wildtypeIDtoCompartment:
				metaboliteIDs.append(
					metaboliteID.upper() + wildtypeIDtoCompartment[metaboliteID]
					)

			else:
				metaboliteIDs.append(
					metaboliteID.upper() + "[c]"
					)

			metaboliteConcentrations.append(concentration)

		# Calculate the following assuming 60 min doubling time

		initWaterMass = sim_data.mass.avgCellWaterMassInit.asNumber(units.g)
		initDryMass = sim_data.mass.avgCellDryMassInit.asNumber(units.g)

		initCellMass = initWaterMass + initDryMass

		initCellVolume = initCellMass / CELL_DENSITY # L

		(massFractions,) = [item for item in raw_data.dryMassComposition if item["doublingTime"] == 60.0 * units.min]

		for entry in raw_data.glycogenFractions:
			metaboliteID = entry["metaboliteId"]

			assert metaboliteID not in metaboliteIDs

			massFrac = entry["massFraction"] * massFractions["glycogenMassFraction"]
			print metaboliteID
			molWeight = sim_data.getter.getMass([metaboliteID])[0].asNumber(units.g / units.mol)

			massInit = massFrac * initDryMass
			molesInit = massInit/molWeight

			concentration = molesInit / initCellVolume

			metaboliteIDs.append(metaboliteID)
			metaboliteConcentrations.append(concentration)

		for entry in raw_data.mureinFractions:
			metaboliteID = entry["metaboliteId"]

			assert metaboliteID not in metaboliteIDs

			massFrac = entry["massFraction"] * massFractions["mureinMassFraction"]
			molWeight = sim_data.getter.getMass([metaboliteID])[0].asNumber(units.g / units.mol)

			massInit = massFrac * initDryMass
			molesInit = massInit/molWeight

			concentration = molesInit / initCellVolume

			metaboliteIDs.append(metaboliteID)
			metaboliteConcentrations.append(concentration)

		for entry in raw_data.LPSFractions:
			metaboliteID = entry["metaboliteId"]

			assert metaboliteID not in metaboliteIDs

			massFrac = entry["massFraction"] * massFractions["lpsMassFraction"]
			molWeight = sim_data.getter.getMass([metaboliteID])[0].asNumber(units.g / units.mol)

			massInit = massFrac * initDryMass
			molesInit = massInit/molWeight

			concentration = molesInit / initCellVolume

			metaboliteIDs.append(metaboliteID)
			metaboliteConcentrations.append(concentration)

		for entry in raw_data.lipidFractions:
			metaboliteID = entry["metaboliteId"]

			assert metaboliteID not in metaboliteIDs

			massFrac = entry["massFraction"] * massFractions["lipidMassFraction"]
			molWeight = sim_data.getter.getMass([metaboliteID])[0].asNumber(units.g / units.mol)

			massInit = massFrac * initDryMass
			molesInit = massInit/molWeight

			concentration = molesInit / initCellVolume

			metaboliteIDs.append(metaboliteID)
			metaboliteConcentrations.append(concentration)

		for entry in raw_data.ionFractions:
			metaboliteID = entry["metaboliteId"]

			assert metaboliteID not in metaboliteIDs

			massFrac = entry["massFraction"] * massFractions["inorganicIonMassFraction"]
			molWeight = sim_data.getter.getMass([metaboliteID])[0].asNumber(units.g / units.mol)

			massInit = massFrac * initDryMass
			molesInit = massInit/molWeight

			concentration = molesInit / initCellVolume

			metaboliteIDs.append(metaboliteID)
			metaboliteConcentrations.append(concentration)

		for entry in raw_data.solubleFractions:
			metaboliteID = entry["metaboliteId"]

			if metaboliteID not in metaboliteIDs:
				massFrac = entry["massFraction"] * massFractions["solublePoolMassFraction"]
				molWeight = sim_data.getter.getMass([metaboliteID])[0].asNumber(units.g / units.mol)

				massInit = massFrac * initDryMass
				molesInit = massInit/molWeight

				concentration = molesInit / initCellVolume

				metaboliteIDs.append(metaboliteID)
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
			metaboliteConcentrations[metaboliteIDs.index("ALA[c]")]
			)

		metaboliteIDs.append("CYS[c]")
		metaboliteConcentrations.append(aaSmallestConc)

		metaboliteIDs.append("SEC[c]")
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

		# H2O: reported water content of E. coli

		h2oMolWeight = sim_data.getter.getMass(["WATER[c]"])[0].asNumber(units.g / units.mol)
		h2oMoles = initWaterMass/h2oMolWeight

		h2oConcentration = h2oMoles / initCellVolume

		metaboliteIDs.append("WATER[c]")
		metaboliteConcentrations.append(h2oConcentration)

		# H: reported pH

		hydrogenConcentration = 10**(-ECOLI_PH)

		metaboliteIDs.append("PROTON[c]")
		metaboliteConcentrations.append(hydrogenConcentration)

		# PPI: multiple sources report 0.5 mM

		# NOTE: Nick says that the physiological levels of PPI are very low - investigate this

		metaboliteIDs.append("PPI[c]")
		metaboliteConcentrations.append(PPI_CONCENTRATION)

		unaccounted = set(wildtypeIDs) - set(metaboliteIDs)

		assert len(unaccounted) == 0

		# Add byproducts with no annotated concentration to force recycling

		metaboliteIDs.append("UMP[c]")
		metaboliteConcentrations.append(2.40e-5)

		# Other quantities to consider:
		# - (d)NTP byproducts not currently included

		self.metabolitePoolIDs = metaboliteIDs
		self.metabolitePoolConcentrations = units.mol/units.L * np.array(metaboliteConcentrations)
		self.cellDensity = units.g/units.L * CELL_DENSITY


	def _buildMetabolism(self, raw_data, sim_data):
		# Build the matrices/vectors for metabolism (FBA)

		# These may be modified/extended later, but should provide the basic
		# data structures

		exchangeReactions = [x["id"] for x in raw_data.reactions if x["id"].startswith("FEIST_EX") or x["id"].startswith("FEIST_DM_")]
		exchangeReactions += ["SELNP_MEDIA_EXCHANGE_HACKED"]

		disabledReactions = ("FEIST_CAT_0", "FEIST_SPODM_0", "FEIST_SPODMpp",
			"FEIST_FHL_0_0", "FEIST_FHL_1_0", "FEIST_ATPM")
		# NOTE: the first five were disabled by Feist, the last is NGAM which model explictly elsewhere

		reactionStoich = {}
		externalExchangeMolecules = []
		reversibleReactions = []
		reactionEnzymes = {}
		reactionRates = {}

		unconstrainedExchangeMolecules = ("CA2[e]", "CL[e]", "CO2[e]", "COBALT2[e]",
			"CU2[e]", "FE2[e]", "FE3[e]", "H[e]", "H2O[e]", "K[e]", "MG2[e]",
			"MN2[e]", "MOBD[e]", "NA1[e]", "NH4[e]", "PI[e]", "SO4[e]", "TUNGS[e]",
			"ZN2[e]", "SELNP[e]",)

		exchangeUnits = units.mmol / units.g / units.h

		constrainedExchangeMolecules = {
			"CBL1[e]": exchangeUnits * 0.01, # TODO: try removing this constraint
			"GLC-D[e]": exchangeUnits * 8,
			"O2[e]": exchangeUnits * 18.5
			}

		catalysisUnits = 1 / units.s

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
			reactionID = reaction["id"]
			stoich = reaction["stoichiometry"]

			if reactionID in disabledReactions:
				continue

			elif reactionID in exchangeReactions:
				if len(stoich) != 1:
					raise Exception("Invalid exchange reaction")

				externalExchangeMolecules.append("{}[{}]".format(
					stoich[0]["molecule"], stoich[0]["location"]
					))

			else:
				if len(stoich) <= 1:
					raise Exception("Invalid biochemical reaction")

				reducedStoich = {
					"{}[{}]".format(
						entry["molecule"], entry["location"]
						): entry["coeff"]
					for entry in stoich
					}

				reactionStoich[reactionID] = reducedStoich

				# Assign reversibilty

				if reaction["dir"] == 0:
					reversibleReactions.append(reactionID)

				# Assign k_cat, if known

				kcat = reaction["kcat"]

				if kcat is not None:
					reactionRates[reactionID] = kcat * catalysisUnits

				# Assign enzyme, if any

				if reactionID in REACTION_ENZYME_ASSOCIATIONS.viewkeys():
					enzymes = REACTION_ENZYME_ASSOCIATIONS[reactionID]

				else:
					enzymes = reaction['catBy']

				reactantLocations = {reactant["location"]
					for reactant in reaction["stoichiometry"]}

				if enzymes is not None and len(enzymes) > 0:
					if len(enzymes) > 1:
						raise Exception("Reaction {} has multiple associated enzymes: {}".format(
							reactionID, enzymes))

					(enzyme,) = enzymes

					validLocations = validEnzymeCompartments[enzyme]

					# TODO: Commenting this so that I don't have to worry about enzymatic complexes. Fix this John.
					# if len(validLocations) == 0:
					# 	raise Exception("Reaction {} uses enzyme {} but this enzyme does not exist.".format(
					# 		reactionID,
					# 		enzyme
					# 		))

					if len(validLocations) == 1:
						(location,) = validLocations

					elif len(reactantLocations) == 1:
						(location,) = reactantLocations

					elif reactantLocations == {"p", "e"}: # if reaction is periplasm <-> extracellular
						(location,) = {"o"} # assume enzyme is in outer membrane

					elif reactantLocations == {"c", "p"}: # if reaction is cytoplasm <-> periplasm
						(location,) = {"i"} # assume enzyme is in inner membrane

					else:
						pass
						# TODO: Commenting this so that I don't have to worry about enzymatic complexes. Fix this John.
						# raise Exception("Reaction {} has multiple associated locations: {}".format(
						# 	reactionID,
						# 	reactantLocations
						# 	))

					# TODO: Commenting this so that I don't have to worry about enzymatic complexes. Fix this John.
					# assert location in validLocations

					enzymeID = "{}[{}]".format(enzyme, location)

					reactionEnzymes[reactionID] = enzymeID

		mws = sim_data.getter.getMass(externalExchangeMolecules)

		exchangeMasses = {moleculeID:mws[index]
			for index, moleculeID in enumerate(externalExchangeMolecules)}

		# Filter out reaction-enzyme associations that lack rates

		reactionEnzymes = {
			reactionID:enzymeID
			for reactionID, enzymeID in reactionEnzymes.viewitems()
			if reactionRates.has_key(reactionID)
			}

		self.reactionStoich = reactionStoich
		self.externalExchangeMolecules = externalExchangeMolecules
		self._exchangeMasses = exchangeMasses
		self.reversibleReactions = reversibleReactions
		self.reactionEnzymes = reactionEnzymes
		self._reactionRates = reactionRates
		self._unconstrainedExchangeMolecules = unconstrainedExchangeMolecules
		self._constrainedExchangeMolecules = constrainedExchangeMolecules


	# def reactionRates(self, timeStep):
	# 	return {
	# 		reactionID:(reactionRate * timeStep).asNumber()
	# 		for reactionID, reactionRate in self._reactionRates.viewitems()
	# 		}


	# def exchangeMasses(self, targetUnits):
	# 	return {
	# 		moleculeID:mass.asNumber(targetUnits)
	# 		for moleculeID, mass in self._exchangeMasses.viewitems()
	# 		}


	def exchangeConstraints(self, exchangeIDs, coefficient, targetUnits):
		externalMoleculeLevels = np.zeros(len(exchangeIDs), np.float64)

		for index, moleculeID in enumerate(exchangeIDs):
			if moleculeID in self._unconstrainedExchangeMolecules:
				externalMoleculeLevels[index] = np.inf

			elif moleculeID in self._constrainedExchangeMolecules.viewkeys():
				externalMoleculeLevels[index] = (
					self._constrainedExchangeMolecules[moleculeID] * coefficient
					).asNumber(targetUnits)

		return externalMoleculeLevels
