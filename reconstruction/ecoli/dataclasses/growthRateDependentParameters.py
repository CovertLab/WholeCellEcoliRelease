"""
SimulationData mass data

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 03/13/2015
"""

from __future__ import division

import numpy as np
from scipy import interpolate
from wholecell.utils import units
import unum

DNA_CRITICAL_MASS = {100: 600, 44: 975} # units of fg
FRACTION_INCREASE_RNAP_PROTEINS = {100: 0, 44: 0.05}

class Mass(object):
	""" Mass """

	def __init__(self, raw_data, sim_data,
				 options):
		self._doubling_time = sim_data.doubling_time

		self._buildConstants(raw_data, sim_data)
		self._buildSubMasses(raw_data, sim_data)
		self._buildCDPeriod(raw_data, sim_data)

		# Alternate masses for paper investigation
		self._alternate_mass_fractions = {
			"protein": options['alternate_mass_fraction_protein'],
			"rna": options['alternate_mass_fraction_rna'],
			"mrna": options['alternate_mass_fraction_mrna']}
		self._buildAlternateSubMasses(raw_data)

		self.avgCellDryMass = self.getAvgCellDryMass(self._doubling_time)

		self._buildDependentConstants()

		self._buildTrnaData(raw_data, sim_data)

	## Setup constants
	def _buildConstants(self, raw_data, sim_data):
		self.__dict__.update(raw_data.mass_parameters)
		self.__dict__.update(raw_data.mass_parameters_alternate)

		self.cellDryMassFraction = 1. - self.cellWaterMassFraction
		self.avgCellToInitialCellConvFactor = np.exp(np.log(2) * self.avgCellCellCycleProgress)

		self._cellDensity = sim_data.constants.cellDensity

		self._glycogenFractions = raw_data.massFractions.glycogenFractions
		self._mureinFractions = raw_data.massFractions.mureinFractions
		self._LPSFractions = raw_data.massFractions.LPSFractions
		self._lipidFractions = raw_data.massFractions.lipidFractions
		self._ionFractions = raw_data.massFractions.ionFractions
		self._solubleFractions = raw_data.massFractions.solubleFractions

		metIds = (
			[x["metaboliteId"] for x in self._glycogenFractions] +
			[x["metaboliteId"] for x in self._mureinFractions] +
			[x["metaboliteId"] for x in self._LPSFractions] +
			[x["metaboliteId"] for x in self._lipidFractions] +
			[x["metaboliteId"] for x in self._ionFractions] +
			[x["metaboliteId"] for x in self._solubleFractions] +
			["WATER[c]"]
			)
		mws = sim_data.getter.getMass(metIds)
		self._mws = dict(zip(metIds, mws))

		self._metTargetIds = [x["Metabolite"] + "[c]" for x in raw_data.metaboliteConcentrations]

	def _buildDependentConstants(self):
		self.avgCellDryMassInit = self.avgCellDryMass / self.avgCellToInitialCellConvFactor
		avgCellWaterMass = (self.avgCellDryMass / self.cellDryMassFraction) * self.cellWaterMassFraction
		self.avgCellWaterMassInit = avgCellWaterMass / self.avgCellToInitialCellConvFactor

	# Setup sub-masses
	def _buildSubMasses(self, raw_data, sim_data):
		self._doubling_time_vector = units.min * np.array([float(x['doublingTime'].asNumber(units.min)) for x in raw_data.dryMassComposition])

		# TODO: Use helper functions written for growthRateDependent parameters to make this better!
		dryMass = np.array([float(x['averageDryMass'].asNumber(units.fg)) for x in raw_data.dryMassComposition])
		self._dryMassParams = interpolate.splrep(self._doubling_time_vector.asNumber(units.min)[::-1], dryMass[::-1])

		self._proteinMassFractionParams = self._getFitParameters(raw_data.dryMassComposition, 'proteinMassFraction')
		self._rnaMassFractionParams = self._getFitParameters(raw_data.dryMassComposition, 'rnaMassFraction')
		self._lipidMassFractionParams = self._getFitParameters(raw_data.dryMassComposition, 'lipidMassFraction')
		self._lpsMassFractionParams = self._getFitParameters(raw_data.dryMassComposition, 'lpsMassFraction')
		self._mureinMassFractionParams = self._getFitParameters(raw_data.dryMassComposition, 'mureinMassFraction')
		self._glycogenMassFractionParams = self._getFitParameters(raw_data.dryMassComposition, 'glycogenMassFraction')
		self._solublePoolMassFractionParams = self._getFitParameters(raw_data.dryMassComposition, 'solublePoolMassFraction')
		self._inorganicIonMassFractionParams = self._getFitParameters(raw_data.dryMassComposition, 'inorganicIonMassFraction')

		self.chromosomeSequenceMass = self._chromosomeSequenceMass(raw_data, sim_data)

	def _buildAlternateSubMasses(self, raw_data):
		'''
		Builds alternate mass fractions according to data described in alternate
		dry mass fraction files:
			dryMassComposition_alternateProtein.tsv
			dryMassComposition_alternateProteinAndRna.tsv
			dryMassComposition_alternateRna.tsv

		In preparation for the interpolation, doubling times that are not
		described in the alternate dry mass fraction files are taken from the
		original dryMassComposition.tsv.
		'''

		base = raw_data.dryMassComposition[:]

		# Alternate protein mass fraction
		alternateProtein = base[:]
		alternateProteinDoublingTime = raw_data.dryMassComposition_alternateProtein[0]["doublingTime"]
		alternateProtein[[x["doublingTime"] for x in base].index(alternateProteinDoublingTime)] = raw_data.dryMassComposition_alternateProtein[0]

		self._proteinMassFractionParams_alternateProtein = self._getFitParameters(alternateProtein, 'proteinMassFraction')
		self._rnaMassFractionParams_alternateProtein = self._getFitParameters(alternateProtein, 'rnaMassFraction')
		self._lipidMassFractionParams_alternateProtein = self._getFitParameters(alternateProtein, 'lipidMassFraction')
		self._lpsMassFractionParams_alternateProtein = self._getFitParameters(alternateProtein, 'lpsMassFraction')
		self._mureinMassFractionParams_alternateProtein = self._getFitParameters(alternateProtein, 'mureinMassFraction')
		self._glycogenMassFractionParams_alternateProtein = self._getFitParameters(alternateProtein, 'glycogenMassFraction')
		self._solublePoolMassFractionParams_alternateProtein = self._getFitParameters(alternateProtein, 'solublePoolMassFraction')
		self._inorganicIonMassFractionParams_alternateProtein = self._getFitParameters(alternateProtein, 'inorganicIonMassFraction')

		# Alternate RNA mass fraction
		alternateRna = base[:]
		alternateRnaDoublingTime = raw_data.dryMassComposition_alternateRna[0]["doublingTime"]
		alternateRna[[x["doublingTime"] for x in base].index(alternateRnaDoublingTime)] = raw_data.dryMassComposition_alternateRna[0]

		self._proteinMassFractionParams_alternateRna = self._getFitParameters(alternateRna, 'proteinMassFraction')
		self._rnaMassFractionParams_alternateRna = self._getFitParameters(alternateRna, 'rnaMassFraction')
		self._lipidMassFractionParams_alternateRna = self._getFitParameters(alternateRna, 'lipidMassFraction')
		self._lpsMassFractionParams_alternateRna = self._getFitParameters(alternateRna, 'lpsMassFraction')
		self._mureinMassFractionParams_alternateRna = self._getFitParameters(alternateRna, 'mureinMassFraction')
		self._glycogenMassFractionParams_alternateRna = self._getFitParameters(alternateRna, 'glycogenMassFraction')
		self._solublePoolMassFractionParams_alternateRna = self._getFitParameters(alternateRna, 'solublePoolMassFraction')
		self._inorganicIonMassFractionParams_alternateRna = self._getFitParameters(alternateRna, 'inorganicIonMassFraction')

		# Alternate protein and RNA mass fractions
		alternateProteinAndRna = base[:]
		alternateProteinAndRnaDoublingTime = raw_data.dryMassComposition_alternateProteinAndRna[0]["doublingTime"]
		alternateProteinAndRna[[x["doublingTime"] for x in base].index(alternateProteinAndRnaDoublingTime)] = raw_data.dryMassComposition_alternateProteinAndRna[0]

		self._proteinMassFractionParams_alternateProteinAndRna = self._getFitParameters(alternateProteinAndRna, 'proteinMassFraction')
		self._rnaMassFractionParams_alternateProteinAndRna = self._getFitParameters(alternateProteinAndRna, 'rnaMassFraction')
		self._lipidMassFractionParams_alternateProteinAndRna = self._getFitParameters(alternateProteinAndRna, 'lipidMassFraction')
		self._lpsMassFractionParams_alternateProteinAndRna = self._getFitParameters(alternateProteinAndRna, 'lpsMassFraction')
		self._mureinMassFractionParams_alternateProteinAndRna = self._getFitParameters(alternateProteinAndRna, 'mureinMassFraction')
		self._glycogenMassFractionParams_alternateProteinAndRna = self._getFitParameters(alternateProteinAndRna, 'glycogenMassFraction')
		self._solublePoolMassFractionParams_alternateProteinAndRna = self._getFitParameters(alternateProteinAndRna, 'solublePoolMassFraction')
		self._inorganicIonMassFractionParams_alternateProteinAndRna = self._getFitParameters(alternateProteinAndRna, 'inorganicIonMassFraction')

	def _getFitParameters(self, dryMassComposition, massFractionName):
		massFraction = np.array([float(x[massFractionName]) for x in dryMassComposition])
		x = self._doubling_time_vector.asNumber(units.min)[::-1]
		y = massFraction[::-1]
		massParams = interpolate.splrep(x, y)
		if np.sum(np.absolute(interpolate.splev(self._doubling_time_vector.asNumber(units.min), massParams) - massFraction)) / massFraction.size > 1.:
			raise Exception("Fitting {} with double exponential, residuals are huge!".format(massFractionName))
		return massParams

	def _buildCDPeriod(self, raw_data, sim_data):
		self.c_period = sim_data.growthRateParameters.c_period
		self.d_period = sim_data.growthRateParameters.d_period

	# Set based on growth rate avgCellDryMass
	def getAvgCellDryMass(self, doubling_time):
		doubling_time = self._clipTau_d(doubling_time)
		avgCellDryMass = units.fg * float(interpolate.splev(doubling_time.asNumber(units.min), self._dryMassParams))
		return avgCellDryMass

	# Set mass fractions based on growth rate
	def getMassFraction(self, doubling_time):
		if type(doubling_time) != unum.Unum:
			raise Exception("Doubling time was not set!")

		D = {}
		dnaMassFraction = self._calculateGrowthRateDependentDnaMass(doubling_time) / self.getAvgCellDryMass(doubling_time)
		dnaMassFraction.normalize()
		dnaMassFraction.checkNoUnit()
		D["dna"] = dnaMassFraction.asNumber()

		doubling_time = self._clipTau_d(doubling_time)

		# Load from a different dry mass composition depending on alternate settings.
		if np.logical_and(self._alternate_mass_fractions["protein"], self._alternate_mass_fractions["rna"]):
			# Alternate protein and RNA mass fractions
			D["protein"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._proteinMassFractionParams_alternateProteinAndRna))
			D["rna"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._rnaMassFractionParams_alternateProteinAndRna))
			D["lipid"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._lipidMassFractionParams_alternateProteinAndRna))
			D["lps"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._lpsMassFractionParams_alternateProteinAndRna))
			D["murein"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._mureinMassFractionParams_alternateProteinAndRna))
			D["glycogen"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._glycogenMassFractionParams_alternateProteinAndRna))
			D["solublePool"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._solublePoolMassFractionParams_alternateProteinAndRna))
			D["inorganicIon"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._inorganicIonMassFractionParams_alternateProteinAndRna))

		elif self._alternate_mass_fractions["protein"]:
			# Alternate protein mass fractions
			D["protein"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._proteinMassFractionParams_alternateProtein))
			D["rna"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._rnaMassFractionParams_alternateProtein))
			D["lipid"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._lipidMassFractionParams_alternateProtein))
			D["lps"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._lpsMassFractionParams_alternateProtein))
			D["murein"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._mureinMassFractionParams_alternateProtein))
			D["glycogen"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._glycogenMassFractionParams_alternateProtein))
			D["solublePool"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._solublePoolMassFractionParams_alternateProtein))
			D["inorganicIon"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._inorganicIonMassFractionParams_alternateProtein))

		elif self._alternate_mass_fractions["rna"]:
			# Alternate RNA mass fractions
			D["protein"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._proteinMassFractionParams_alternateRna))
			D["rna"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._rnaMassFractionParams_alternateRna))
			D["lipid"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._lipidMassFractionParams_alternateRna))
			D["lps"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._lpsMassFractionParams_alternateRna))
			D["murein"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._mureinMassFractionParams_alternateRna))
			D["glycogen"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._glycogenMassFractionParams_alternateRna))
			D["solublePool"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._solublePoolMassFractionParams_alternateRna))
			D["inorganicIon"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._inorganicIonMassFractionParams_alternateRna))

		else:
			# Original mass fractions
			D["protein"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._proteinMassFractionParams))
			D["rna"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._rnaMassFractionParams))
			D["lipid"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._lipidMassFractionParams))
			D["lps"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._lpsMassFractionParams))
			D["murein"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._mureinMassFractionParams))
			D["glycogen"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._glycogenMassFractionParams))
			D["solublePool"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._solublePoolMassFractionParams))
			D["inorganicIon"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._inorganicIonMassFractionParams))

		total = np.sum([y for x,y in D.iteritems()])
		for key, value in D.iteritems():
			if key != 'dna':
				D[key] = value / total
		assert np.absolute(np.sum([x for x in D.itervalues()]) - 1.) < 1e-3
		return D


	def getFractionMass(self, doubling_time):
		D = {}
		massFraction = self.getMassFraction(doubling_time)
		for key, value in massFraction.iteritems():
			D[key + "Mass"] = value * self.getAvgCellDryMass(doubling_time)

		# Load from a different RNA mass fraction depending on alternate settings.
		if self._alternate_mass_fractions["mrna"]:
			D["rRna23SMass"] = D['rnaMass'] * self._rrna23s_mass_sub_fraction_alternate
			D["rRna16SMass"] = D['rnaMass'] * self._rrna16s_mass_sub_fraction_alternate
			D["rRna5SMass"] = D['rnaMass'] * self._rrna5s_mass_sub_fraction_alternate
			D["tRnaMass"] = D['rnaMass'] * self._trna_mass_sub_fraction_alternate
			D["mRnaMass"] = D['rnaMass'] * self._mrna_mass_sub_fraction_alternate

		else:
			D["rRna23SMass"] = D['rnaMass'] * self._rrna23s_mass_sub_fraction
			D["rRna16SMass"] = D['rnaMass'] * self._rrna16s_mass_sub_fraction
			D["rRna5SMass"] = D['rnaMass'] * self._rrna5s_mass_sub_fraction
			D["tRnaMass"] = D['rnaMass'] * self._trna_mass_sub_fraction
			D["mRnaMass"] = D['rnaMass'] * self._mrna_mass_sub_fraction

		return D

	def getBiomassAsConcentrations(self, doubling_time):

		avgCellDryMassInit = self.getAvgCellDryMass(doubling_time) / self.avgCellToInitialCellConvFactor
		avgCellWaterMassInit = (avgCellDryMassInit / self.cellDryMassFraction) * self.cellWaterMassFraction

		initWaterMass = avgCellWaterMassInit.asNumber(units.g)
		initDryMass = avgCellDryMassInit.asNumber(units.g)

		initCellMass = initWaterMass + initDryMass

		initCellVolume = initCellMass / self._cellDensity.asNumber(units.g / units.L) # L

		massFraction = self.getMassFraction(doubling_time)

		metaboliteIDs = []
		metaboliteConcentrations = []

		for entry in self._glycogenFractions:
			metaboliteID = entry["metaboliteId"]

			assert metaboliteID not in metaboliteIDs + self._metTargetIds

			massFrac = entry["massFraction"] * massFraction["glycogen"]
			molWeight = self._mws[metaboliteID].asNumber(units.g / units.mol)

			massInit = massFrac * initDryMass
			molesInit = massInit / molWeight

			concentration = molesInit / initCellVolume

			metaboliteIDs.append(metaboliteID)
			metaboliteConcentrations.append(concentration)

		for entry in self._mureinFractions:
			metaboliteID = entry["metaboliteId"]

			assert metaboliteID not in metaboliteIDs + self._metTargetIds

			massFrac = entry["massFraction"] * massFraction["murein"]
			molWeight = self._mws[metaboliteID].asNumber(units.g / units.mol)

			massInit = massFrac * initDryMass
			molesInit = massInit / molWeight

			concentration = molesInit / initCellVolume

			metaboliteIDs.append(metaboliteID)
			metaboliteConcentrations.append(concentration)

		for entry in self._LPSFractions:
			metaboliteID = entry["metaboliteId"]

			assert metaboliteID not in metaboliteIDs + self._metTargetIds

			massFrac = entry["massFraction"] * massFraction["lps"]
			molWeight = self._mws[metaboliteID].asNumber(units.g / units.mol)

			massInit = massFrac * initDryMass
			molesInit = massInit / molWeight

			concentration = molesInit / initCellVolume

			metaboliteIDs.append(metaboliteID)
			metaboliteConcentrations.append(concentration)

		for entry in self._lipidFractions:
			metaboliteID = entry["metaboliteId"]

			assert metaboliteID not in metaboliteIDs + self._metTargetIds

			massFrac = entry["massFraction"] * massFraction["lipid"]
			molWeight = self._mws[metaboliteID].asNumber(units.g / units.mol)

			massInit = massFrac * initDryMass
			molesInit = massInit / molWeight

			concentration = molesInit / initCellVolume

			metaboliteIDs.append(metaboliteID)
			metaboliteConcentrations.append(concentration)

		for entry in self._ionFractions:
			metaboliteID = entry["metaboliteId"]

			assert metaboliteID not in metaboliteIDs + self._metTargetIds

			massFrac = entry["massFraction"] * massFraction["inorganicIon"]
			molWeight = self._mws[metaboliteID].asNumber(units.g / units.mol)

			massInit = massFrac * initDryMass
			molesInit = massInit / molWeight

			concentration = molesInit / initCellVolume

			metaboliteIDs.append(metaboliteID)
			metaboliteConcentrations.append(concentration)

		for entry in self._solubleFractions:
			metaboliteID = entry["metaboliteId"]

			if metaboliteID not in self._metTargetIds:
				massFrac = entry["massFraction"] * massFraction["solublePool"]
				molWeight = self._mws[metaboliteID].asNumber(units.g / units.mol)

				massInit = massFrac * initDryMass
				molesInit = massInit / molWeight

				concentration = molesInit / initCellVolume

				metaboliteIDs.append(metaboliteID)
				metaboliteConcentrations.append(concentration)


		# H2O: reported water content of E. coli

		h2oMolWeight = self._mws["WATER[c]"].asNumber(units.g / units.mol)
		h2oMoles = initWaterMass / h2oMolWeight

		h2oConcentration = h2oMoles / initCellVolume

		metaboliteIDs.append("WATER[c]")
		metaboliteConcentrations.append(h2oConcentration)


		metaboliteConcentrations = (units.mol / units.L) * np.array(metaboliteConcentrations)

		return dict(zip(metaboliteIDs, metaboliteConcentrations))

	def _calculateGrowthRateDependentDnaMass(self, doubling_time):
		C_PERIOD = self.c_period
		D_PERIOD = self.d_period
		CD_PERIOD = C_PERIOD + D_PERIOD

		if doubling_time < D_PERIOD:
			raise Exception, "Can't have doubling time shorter than cytokinesis time!"

		doubling_time_unit = units.getUnit(doubling_time)

		# TODO: If you really care, this should be a loop.
		# It is optimized to run quickly over the range of T_d
		# and C and D periods that we have.
		return self.chromosomeSequenceMass * (1 +
			1 * (np.maximum(0. * doubling_time_unit, CD_PERIOD - doubling_time) / C_PERIOD) +
			2 * (np.maximum(0. * doubling_time_unit, CD_PERIOD - 2 * doubling_time) / C_PERIOD) +
			4 * (np.maximum(0. * doubling_time_unit, CD_PERIOD - 4 * doubling_time) / C_PERIOD)
			)

	def _chromosomeSequenceMass(self, raw_data, sim_data):
		dntCounts = np.array([
			raw_data.genome_sequence.count('A') + raw_data.genome_sequence.count('T'),
			raw_data.genome_sequence.count('C') + raw_data.genome_sequence.count('G'),
			raw_data.genome_sequence.count('G') + raw_data.genome_sequence.count('C'),
			raw_data.genome_sequence.count('T') + raw_data.genome_sequence.count('A')
		])

		dntMasses = (sim_data.getter.getMass(sim_data.moleculeGroups.polymerizedDNT_IDs) / raw_data.constants['nAvogadro']).asUnit(units.g)
		chromMass = units.dot(dntCounts, dntMasses)
		return chromMass

	def _clipTau_d(self, doubling_time):
		# Clip values to be in the range that we have data for
		if hasattr(doubling_time, "dtype"):
			doubling_time[doubling_time > max(self._doubling_time_vector)] = max(self._doubling_time_vector)
			doubling_time[doubling_time < min(self._doubling_time_vector)] = min(self._doubling_time_vector)
		else:
			if doubling_time > max(self._doubling_time_vector):
				doubling_time = max(self._doubling_time_vector)
			elif doubling_time < min(self._doubling_time_vector):
				doubling_time = min(self._doubling_time_vector)
		return doubling_time

	def _buildTrnaData(self, raw_data, sim_data):
		growth_rate_unit = units.getUnit(raw_data.trnaData.trna_growth_rates[0]['growth rate'])

		self._trna_growth_rates = growth_rate_unit * np.array([x['growth rate'].asNumber() for x in raw_data.trnaData.trna_growth_rates])

		trna_ratio_to_16SrRNA_by_growth_rate = []
		for gr in self._trna_growth_rates: # This is a little crazy...
			trna_ratio_to_16SrRNA_by_growth_rate.append([x['ratio to 16SrRNA'] for x in getattr(raw_data.trnaData, "trna_ratio_to_16SrRNA_" + str(gr.asNumber()).replace('.','p'))])
		self._trna_ratio_to_16SrRNA_by_growth_rate = np.array(trna_ratio_to_16SrRNA_by_growth_rate)

		self._trna_ids = [x['rna id'] for x in raw_data.trnaData.trna_ratio_to_16SrRNA_0p4]

	def getTrnaDistribution(self):
		return self._getTrnaAbundanceAtGrowthRate(self._doubling_time)

	def _getTrnaAbundanceAtGrowthRate(self, doubling_time):
		assert type(doubling_time) == unum.Unum
		assert type(doubling_time.asNumber()) == float
		growth_rate = 1 / doubling_time
		growth_rate = growth_rate.asNumber(1/units.h)

		from scipy.interpolate import interp1d
		trna_abundance_interpolation_functions = [interp1d(self._trna_growth_rates.asNumber(1/units.h), self._trna_ratio_to_16SrRNA_by_growth_rate[:,i]) for i in range(self._trna_ratio_to_16SrRNA_by_growth_rate.shape[1])]

		abundance = np.zeros(len(self._trna_ids), dtype = [('id','a50'),('molar_ratio_to_16SrRNA', np.float64)])
		abundance['id'] = self._trna_ids
		abundance['molar_ratio_to_16SrRNA'] = [x(growth_rate) for x in trna_abundance_interpolation_functions]
		return abundance

class GrowthRateParameters(object):
	"""
	GrowthRateParameters
	"""

	def __init__(self,
				 raw_data,
				 sim_data,
				 options):
		self._doubling_time = sim_data.doubling_time

		if options['alternate_ribosome_activity']:
			_loadTableIntoObjectGivenDoublingTime(self, raw_data.growthRateDependentParameters_alternateRibosomeActivity)
			self.ribosomeElongationRateParams = _getFitParameters(raw_data.growthRateDependentParameters_alternateRibosomeActivity, "ribosomeElongationRate")
			self.rnaPolymeraseElongationRateParams = _getFitParameters(raw_data.growthRateDependentParameters_alternateRibosomeActivity, "rnaPolymeraseElongationRate")
			self.fractionActiveRnapParams = _getFitParameters(raw_data.growthRateDependentParameters_alternateRibosomeActivity, "fractionActiveRnap")
			self.fractionActiveRibosomeParams = _getFitParameters(raw_data.growthRateDependentParameters_alternateRibosomeActivity, "fractionActiveRibosome")
		else:
			_loadTableIntoObjectGivenDoublingTime(self, raw_data.growthRateDependentParameters)
			self.ribosomeElongationRateParams = _getFitParameters(raw_data.growthRateDependentParameters, "ribosomeElongationRate")
			self.rnaPolymeraseElongationRateParams = _getFitParameters(raw_data.growthRateDependentParameters, "rnaPolymeraseElongationRate")
			self.fractionActiveRnapParams = _getFitParameters(raw_data.growthRateDependentParameters, "fractionActiveRnap")
			self.fractionActiveRibosomeParams = _getFitParameters(raw_data.growthRateDependentParameters, "fractionActiveRibosome")

			# Replace alternate values
			self.fractionActiveRnap = self.getFractionActiveRnap(self._doubling_time)

		self.c_period = units.min * 40.
		self.d_period = units.min * 20.
		self.dnaPolymeraseElongationRate = units.nt / units.s * 967.

	def getRibosomeElongationRate(self, doubling_time):
		return _useFitParameters(doubling_time, **self.ribosomeElongationRateParams)

	def getRnapElongationRate(self, doubling_time):
		return _useFitParameters(doubling_time, **self.rnaPolymeraseElongationRateParams)

	def getFractionActiveRnap(self, doubling_time):
		return _useFitParameters(doubling_time, **self.fractionActiveRnapParams)

	def getFractionActiveRibosome(self, doubling_time):
		return _useFitParameters(doubling_time, **self.fractionActiveRibosomeParams)

	def getDnaCriticalMass(self, doubling_time):
		return DNA_CRITICAL_MASS.get(doubling_time.asNumber(units.min), DNA_CRITICAL_MASS[44]) * units.fg

	def getFractionIncreaseRnapProteins(self, doubling_time):
		return FRACTION_INCREASE_RNAP_PROTEINS.get(doubling_time.asNumber(units.min), FRACTION_INCREASE_RNAP_PROTEINS[44])

def _getFitParameters(list_of_dicts, key):
	# Load rows of data
	x = _loadRow('doublingTime', list_of_dicts)
	y = _loadRow(key, list_of_dicts)

	# Save and strip units
	y_units = 1
	x_units = 1
	if units.hasUnit(y):
		y_units = units.getUnit(y)
		y = y.asNumber(y_units)
	if units.hasUnit(x):
		x_units = units.getUnit(x)
		x = x.asNumber(x_units)

	# Sort data for spine fitting (must be ascending order)
	idx_order = x.argsort()
	x = x[idx_order]
	y = y[idx_order]

	# Generate fit
	parameters = interpolate.splrep(x, y)
	if np.sum(np.absolute(interpolate.splev(x, parameters) - y)) / y.size > 1.:
		raise Exception("Fitting {} with 3d spline, residuals are huge!".format(key))

	return {'parameters' : parameters, 'x_units' : x_units, 'y_units' : y_units, 'dtype' : y.dtype}

def _useFitParameters(x_new, parameters, x_units, y_units, dtype):
	# Convert to same unit base
	if units.hasUnit(x_units):
		x_new = x_new.asNumber(x_units)
	elif units.hasUnit(x_new):
		raise Exception("New x value has units but fit does not!")

	# Calculate new interpolated y value
	y_new = interpolate.splev(x_new, parameters)

	# If value should be an integer (i.e. an elongation rate)
	# round to the nearest integer
	if dtype == np.int:
		y_new = int(np.round(y_new))

	return y_units * y_new

def _loadRow(key, list_of_dicts):
	if units.hasUnit(list_of_dicts[0][key]):
		row_units = units.getUnit(list_of_dicts[0][key])
		return row_units * np.array([x[key].asNumber(row_units) for x in list_of_dicts])
	else:
		return np.array([x[key] for x in list_of_dicts])

def _loadTableIntoObjectGivenDoublingTime(obj, list_of_dicts):
	table_keys = list_of_dicts[0].keys()

	if 'doublingTime' not in table_keys:
		raise Exception, 'This data has no doubling time column but it is supposed to be growth rate dependent!'
	else:
		table_keys.pop(table_keys.index('doublingTime'))

	for key in table_keys:
		fitParameters = _getFitParameters(list_of_dicts, key)
		attrValue = _useFitParameters(obj._doubling_time, **fitParameters)
		setattr(obj, key, attrValue)
