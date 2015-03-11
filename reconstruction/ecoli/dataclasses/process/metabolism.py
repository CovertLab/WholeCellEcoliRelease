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

CELL_DENSITY = 1.1e3 # g/L
# from Baldwin WW, Myer R, Powell N, Anderson E, Koch AL. Buoyant
# density of Escherichia coli is determined solely by the osmolarity
# of the culture medium. Arch Microbiol. 1995 Aug164(2):155-7 p.156
# fig.3 & fig.2 (retrieved from Bionumbers)

METABOLITE_CONCENTRATIONS = { # mol / L # TODO: move to SQL # TODO: Add units!
	"glu-L": 9.60e-2,
	"gthrd": 1.70e-2,
	"fdp": 1.50e-2,
	"atp": 9.60e-3,
	"u3aga": 9.20e-3,
	"utp": 8.30e-3,
	"gtp": 4.90e-3,
	"dttp": 4.60e-3,
	"asp-L": 4.20e-3,
	"val-L": 4.00e-3,
	"6pgc": 3.80e-3,
	"gln-L": 3.80e-3,
	"ctp": 2.70e-3,
	"ala-L": 2.60e-3,
	"nad": 2.60e-3,
	"udpg": 2.50e-3,
	"uri": 2.10e-3,
	"cit": 2.00e-3,
	"udp": 1.80e-3,
	"mal-L": 1.70e-3,
	"3pg": 1.50e-3,
	"citr-L": 1.40e-3,
	"coa": 1.40e-3,
	"glyc-R": 1.40e-3,
	"gam6p": 1.20e-3,
	"actp": 1.10e-3,
	"6pgl": 1.00e-3,
	"gdp": 6.80e-4,
	"accoa": 6.10e-4,
	"cbasp": 5.90e-4,
	"arg-L": 5.70e-4,
	"succ": 5.70e-4,
	"udpglcur": 5.70e-4,
	"adp": 5.60e-4,
	"asn-L": 5.10e-4,
	"akg": 4.40e-4,
	"lys-L": 4.10e-4,
	"pro-L": 3.90e-4,
	"dtdp": 3.80e-4,
	"dhap": 3.70e-4,
	"hcys-L": 3.70e-4,
	"cmp": 3.60e-4,
	"amp": 2.80e-4,
	"succoa": 2.30e-4,
	"gua": 1.90e-4,
	"pep": 1.80e-4,
	"amet": 1.80e-4,
	"thr-L": 1.80e-4,
	"fad": 1.70e-4,
	"met-L": 1.50e-4,
	"23dhb": 1.40e-4,
	"fum": 1.20e-4,
	"nadph": 1.20e-4,
	"phpyr": 9.00e-5,
	"nadh": 8.30e-5,
	"acgam1p": 8.20e-5,
	"his-L": 6.80e-5,
	"ser-L": 6.80e-5,
	"4hbz": 5.20e-5,
	"dgmp": 5.10e-5,
	"glyc3p": 4.90e-5,
	"acorn": 4.30e-5,
	"glcn": 4.20e-5,
	# "23camp": 3.50e-5, # can't be formed by the reaction network
	"dctp": 3.50e-5,
	"malcoa": 3.50e-5,
	"tyr-L": 2.90e-5,
	"gmp": 2.40e-5,
	"aacoa": 2.20e-5,
	"ribflv": 1.90e-5,
	"phe-L": 1.80e-5,
	"acon-C": 1.60e-5,
	"datp": 1.60e-5,
	"csn": 1.40e-5,
	"skm": 1.40e-5,
	"histd": 1.30e-5,
	"dhor-S": 1.20e-5,
	"quln": 1.20e-5,
	"trp-L": 1.20e-5,
	"orn": 1.00e-5,
	"damp": 8.80e-6,
	"aps": 6.60e-6,
	# "inost": 5.70e-6, # can't be formed by the reaction network
	"ppcoa": 5.30e-6,
	"adpglc": 4.30e-6,
	"anth": 3.50e-6,
	"dad-2": 2.80e-6,
	"cytd": 2.60e-6,
	"nadp": 2.10e-6,
	"gsn": 1.60e-6,
	"ade": 1.50e-6,
	"dgsn": 5.20e-7,
	"adn": 1.30e-7,
	}

ILE_LEU_CONCENTRATION = 3.0e-4 # mmol/L
ILE_FRACTION = 0.360 # the fraction of iso/leucine that is isoleucine; computed from our monomer data
ECOLI_PH = 7.2

PPI_CONCENTRATION = 0.5e-3 # M, multiple sources

class Metabolism(object):
	""" Metabolism """

	def __init__(self, raw_data, sim_data):
		self._buildBiomass(raw_data, sim_data)


	def _buildBiomass(self, raw_data, sim_data):
		coreBiomassIds = []
		coreBiomassConcentration = []
		for metabolite in raw_data.metabolites:
			for location in metabolite['core_location']:
				coreBiomassIds.append('{}[{}]'.format(metabolite['id'], location))
			for concentration in metabolite['core_concentration']:
				coreBiomassConcentration.append(concentration)

		wildtypeBiomassIds = []
		wildtypeBiomassConcentration = []
		for metabolite in raw_data.metabolites:
			for location in metabolite['wildtype_location']:
				wildtypeBiomassIds.append('{}[{}]'.format(metabolite['id'], location))
			for concentration in metabolite['wildtype_concentration']:
				wildtypeBiomassConcentration.append(concentration)

		coreBiomassData = np.zeros(len(coreBiomassIds), dtype = [('id', 'a50'), ('biomassFlux', 'f8')])
		wildtypeBiomassData = np.zeros(len(wildtypeBiomassIds), dtype = [('id', 'a50'), ('biomassFlux',		'f8')])

		coreBiomassData['id'] = coreBiomassIds
		coreBiomassData['biomassFlux'] = coreBiomassConcentration
		wildtypeBiomassData['id'] = wildtypeBiomassIds
		wildtypeBiomassData['biomassFlux'] = wildtypeBiomassConcentration

		field_units = {'id' : None,
				'biomassFlux' : units.mmol / units.g}
		self.coreBiomass 		= UnitStructArray(coreBiomassData, field_units)
		self.wildtypeBiomass 	= UnitStructArray(wildtypeBiomassData, field_units)

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

		wildtypeIDs = wildtypeBiomassData["id"].tolist()

		wildtypeIDtoCompartment = {
			wildtypeID[:-3] : wildtypeID[-3:]
			for wildtypeID in wildtypeIDs
			}

		for metaboliteID, concentration in METABOLITE_CONCENTRATIONS.viewitems():
			if metaboliteID in wildtypeIDtoCompartment:
				metaboliteIDs.append(
					metaboliteID.upper() + wildtypeIDtoCompartment[metaboliteID]
					)

			elif metaboliteID == "23CAMP":
				metaboliteIDs.append(
					metaboliteID.upper() + "[p]"
					)

			else:
				metaboliteIDs.append(
					metaboliteID.upper() + "[c]"
					)

			metaboliteConcentrations.append(concentration)

		# Calculate the following assuming 60 min doubling time

		initWaterMass = raw_data.parameters["avgCellWaterMassInit"].asNumber(units.g)
		initDryMass = raw_data.parameters["avgCellDryMassInit"].asNumber(units.g)

		initCellMass = initWaterMass + initDryMass

		initCellVolume = initCellMass / CELL_DENSITY # L

		(massFractions,) = [item for item in raw_data.dryMassComposition if item["doublingTime"] == 60.0]

		for entry in raw_data.glycogenFractions:
			metaboliteID = entry["metaboliteId"]

			assert metaboliteID not in metaboliteIDs

			massFrac = entry["massFraction"] * massFractions["glycogenMassFraction"]
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

		metaboliteIDs.append("ILE-L[c]")
		metaboliteConcentrations.append(ileRelative * ILE_LEU_CONCENTRATION)

		metaboliteIDs.append("LEU-L[c]")
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
			metaboliteConcentrations[metaboliteIDs.index("ALA-L[c]")]
			)

		metaboliteIDs.append("CYS-L[c]")
		metaboliteConcentrations.append(aaSmallestConc)

		metaboliteIDs.append("SEC-L[c]")
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

		h2oMolWeight = sim_data.getter.getMass(["H2O[c]"])[0].asNumber(units.g / units.mol)
		h2oMoles = initWaterMass/h2oMolWeight

		h2oConcentration = h2oMoles / initCellVolume

		metaboliteIDs.append("H2O[c]")
		metaboliteConcentrations.append(h2oConcentration)

		# H: reported pH

		hydrogenConcentration = 10**(-ECOLI_PH)

		metaboliteIDs.append("H[c]")
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

