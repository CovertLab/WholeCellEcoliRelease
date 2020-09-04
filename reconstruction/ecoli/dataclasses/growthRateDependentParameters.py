"""
SimulationData mass data

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 03/13/2015
"""

from __future__ import absolute_import, division, print_function

from typing import Tuple

import numpy as np
from scipy import interpolate, stats
import unum

from wholecell.utils import units
from six.moves import range
import six
from six.moves import zip


NORMAL_CRITICAL_MASS = 975 * units.fg
SLOW_GROWTH_FACTOR = 1.2  # adjustment for smaller cells


class Mass(object):
	""" Mass """

	def __init__(self, raw_data, sim_data):
		self._doubling_time = sim_data.doubling_time

		self._buildConstants(raw_data, sim_data)
		self._buildSubMasses(raw_data, sim_data)
		self._buildCDPeriod(raw_data, sim_data)

		self.avgCellDryMass = self.getAvgCellDryMass(self._doubling_time)
		self.massFraction = self.getMassFraction(self._doubling_time)
		self.avgCellSubMass = self.getFractionMass(self._doubling_time)

		self._buildDependentConstants()

		self._buildTrnaData(raw_data, sim_data)

	## Setup constants
	def _buildConstants(self, raw_data, sim_data):
		mass_parameters = raw_data.mass_parameters
		self.__dict__.update(mass_parameters)

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

		dryMass = np.array([
			float(x['averageDryMass'].asNumber(units.fg))
			for x in raw_data.dryMassComposition
			])
		self._dryMassParams = linear_regression(
			self._doubling_time_vector.asNumber(units.min), 1. / dryMass)

		self._proteinMassFractionParams = self._getFitParameters(raw_data.dryMassComposition, 'proteinMassFraction')
		self._rnaMassFractionParams = self._getFitParameters(raw_data.dryMassComposition, 'rnaMassFraction')
		self._lipidMassFractionParams = self._getFitParameters(raw_data.dryMassComposition, 'lipidMassFraction')
		self._lpsMassFractionParams = self._getFitParameters(raw_data.dryMassComposition, 'lpsMassFraction')
		self._mureinMassFractionParams = self._getFitParameters(raw_data.dryMassComposition, 'mureinMassFraction')
		self._glycogenMassFractionParams = self._getFitParameters(raw_data.dryMassComposition, 'glycogenMassFraction')
		self._solublePoolMassFractionParams = self._getFitParameters(raw_data.dryMassComposition, 'solublePoolMassFraction')
		self._inorganicIonMassFractionParams = self._getFitParameters(raw_data.dryMassComposition, 'inorganicIonMassFraction')

		self.chromosomeSequenceMass = (
			sim_data.getter.getMass([sim_data.moleculeIds.full_chromosome])[0]
				/ sim_data.constants.nAvogadro
			).asUnit(units.g)

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
		# type: (units.Unum) -> units.Unum
		"""
		Gets the dry mass for an average cell at the given doubling time.

		Args:
			doubling_time (float, time units): expected doubling time

		Returns:
			average cell dry mass (float, mass units)
		"""

		doubling_time = doubling_time.asNumber(units.min)
		inverse_mass = self._dryMassParams[0] * doubling_time + self._dryMassParams[1]
		if inverse_mass < 0:
			raise ValueError('Doubling time ({} min) is too short, could not get mass.'
				.format(doubling_time))
		return units.fg / inverse_mass

	def get_dna_critical_mass(self, doubling_time):
		# type: (units.Unum) -> units.Unum
		"""
		Returns the critical mass for replication initiation.  Faster growing
		cells maintain a consistent initiation mass but slower growing cells
		are smaller and will never reach this mass so it needs to be adjusted
		lower for them.

		Args:
			doubling_time (float with time units): expected doubling time of cell

		Returns:
			critical_mass (float with mass units): critical mass for DNA
				replication initiation
		"""

		mass = self.getAvgCellDryMass(doubling_time) / self.cellDryMassFraction
		critical_mass = min(mass * SLOW_GROWTH_FACTOR, NORMAL_CRITICAL_MASS)
		return critical_mass

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
		D["protein"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._proteinMassFractionParams))
		D["rna"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._rnaMassFractionParams))
		D["lipid"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._lipidMassFractionParams))
		D["lps"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._lpsMassFractionParams))
		D["murein"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._mureinMassFractionParams))
		D["glycogen"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._glycogenMassFractionParams))
		D["solublePool"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._solublePoolMassFractionParams))
		D["inorganicIon"] = float(interpolate.splev(doubling_time.asNumber(units.min), self._inorganicIonMassFractionParams))

		total = np.sum([y for y in six.viewvalues(D)])
		for key, value in six.viewitems(D):
			if key != 'dna':
				D[key] = value / total
		assert np.absolute(np.sum([x for x in six.viewvalues(D)]) - 1.) < 1e-3
		return D


	def getFractionMass(self, doubling_time):
		D = {}
		massFraction = self.getMassFraction(doubling_time)
		for key, value in six.viewitems(massFraction):
			D[key + "Mass"] = value * self.getAvgCellDryMass(doubling_time)

		return D

	def get_basal_rna_fractions(self):
		"""
		Measured RNA subgroup mass fractions. Fractions should change in other
		conditions with growth rate (see transcription.get_rna_fractions()).
		"""

		return {
			'23S': self._rrna23s_mass_sub_fraction,
			'16S': self._rrna16s_mass_sub_fraction,
			'5S': self._rrna5s_mass_sub_fraction,
			'trna': self._trna_mass_sub_fraction,
			'mrna': self._mrna_mass_sub_fraction,
			}

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
			raise Exception(
				"Can't have doubling time shorter than cytokinesis time!")

		doubling_time_unit = units.getUnit(doubling_time)

		# TODO: If you really care, this should be a loop.
		# It is optimized to run quickly over the range of T_d
		# and C and D periods that we have.
		return self.chromosomeSequenceMass * (1 +
			1 * (np.maximum(0. * doubling_time_unit, CD_PERIOD - doubling_time) / C_PERIOD) +
			2 * (np.maximum(0. * doubling_time_unit, CD_PERIOD - 2 * doubling_time) / C_PERIOD) +
			4 * (np.maximum(0. * doubling_time_unit, CD_PERIOD - 4 * doubling_time) / C_PERIOD)
			)

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

	def getTrnaDistribution(self, doubling_time):
		assert type(doubling_time) == unum.Unum
		assert type(doubling_time.asNumber()) == float
		growth_rate = 1 / doubling_time
		growth_rate = growth_rate.asNumber(1/units.h)

		trna_abundance_interpolation_functions = [
			interpolate.interp1d(self._trna_growth_rates.asNumber(1/units.h), self._trna_ratio_to_16SrRNA_by_growth_rate[:,i])
			for i in range(self._trna_ratio_to_16SrRNA_by_growth_rate.shape[1])]

		id_length = max(len(id_) for id_ in self._trna_ids)
		abundance = np.zeros(
			len(self._trna_ids),
			dtype=[
				('id','U{}'.format(id_length)),
				('molar_ratio_to_16SrRNA', np.float64),
			])
		abundance['id'] = self._trna_ids
		abundance['molar_ratio_to_16SrRNA'] = [x(growth_rate) for x in trna_abundance_interpolation_functions]
		return abundance

class GrowthRateParameters(object):
	"""
	GrowthRateParameters
	"""

	def __init__(self, raw_data, sim_data):
		self._doubling_time = sim_data.doubling_time
		_loadTableIntoObjectGivenDoublingTime(self, raw_data.growthRateDependentParameters)
		self.ribosomeElongationRateParams = _getFitParameters(raw_data.growthRateDependentParameters, "ribosomeElongationRate")
		self.rnaPolymeraseElongationRateParams = _getFitParameters(raw_data.growthRateDependentParameters, "rnaPolymeraseElongationRate")
		self.fractionActiveRnapParams = _getFitParameters(raw_data.growthRateDependentParameters, "fractionActiveRnap")
		self.fractionActiveRibosomeParams = _getFitParameters(raw_data.growthRateDependentParameters, "fractionActiveRibosome")
		self.ppGppConcentration = _getFitParameters(raw_data.growthRateDependentParameters, "ppGpp_conc")

		self._per_dry_mass_to_per_volume = sim_data.constants.cellDensity * (1. - raw_data.mass_parameters['cellWaterMassFraction'])
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

	def getppGppConc(self, doubling_time):
		return _useFitParameters(doubling_time, **self.ppGppConcentration) * self._per_dry_mass_to_per_volume

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
	cs = interpolate.CubicSpline(x, y, bc_type='natural')
	if np.sum(np.absolute(cs(x) - y)) / y.size > 1.:
		raise Exception("Fitting {} with 3d spline, residuals are huge!".format(key))

	return {'function': cs, 'x_units': x_units, 'y_units': y_units, 'dtype': y.dtype}

def _useFitParameters(x_new, function, x_units, y_units, dtype):
	# Convert to same unit base
	if units.hasUnit(x_units):
		x_new = x_new.asNumber(x_units)
	elif units.hasUnit(x_new):
		raise Exception("New x value has units but fit does not!")

	# Calculate new interpolated y value
	y_new = function(x_new)

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
	table_keys = list(list_of_dicts[0].keys())

	if 'doublingTime' not in table_keys:
		raise Exception(
			'This data has no doubling time column but it is supposed to be growth rate dependent!')
	else:
		table_keys.pop(table_keys.index('doublingTime'))

	for key in table_keys:
		fitParameters = _getFitParameters(list_of_dicts, key)
		attrValue = _useFitParameters(obj._doubling_time, **fitParameters)
		setattr(obj, key, attrValue)

def linear_regression(x, y, r_tol=0.999, p_tol=1e-5):
	# type: (np.ndarray, np.ndarray, float, float) -> Tuple[float, float]
	"""
	Perform linear regression on a data set and check that statistics are
	within expected values to confirm a good linear fit.

	Args:
		x (float): x values for regression
		y (float): y values for regression
		r_tol: lower limit for r statistic
		p_tol: upper limit for p statistic

	Returns:
		slope: linear fit slope
		intercept: linear fit intercept
	"""

	result = stats.linregress(x, y)

	if result.rvalue < r_tol:
		raise ValueError('Could not fit linear regression with high enough r'
			' value: {} < {}'.format(result.rvalue, r_tol))
	if result.pvalue > p_tol:
		raise ValueError('Could not fit linear regression with low enough p'
			' value: {} > {}'.format(result.pvalue, p_tol))

	return result.slope, result.intercept
