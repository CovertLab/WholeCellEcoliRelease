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

class Mass(object):
	""" Mass """

	def __init__(self, raw_data, sim_data):
		self._doubling_time = sim_data.doubling_time

		self._buildConstants(raw_data, sim_data)
		self._buildSubMasses(raw_data, sim_data)
		self._buildCDPeriod(raw_data, sim_data)

		self.avgCellDryMass = self._setAvgCellDryMass()
		self.avgCell60MinDoublingTimeTotalMassInit = 813.248540675 * units.fg
		self.massFraction = self._setMassFraction()
		self.subMass = self._setSubMass()

		self._buildDependentConstants()

		self._buildTrnaData(raw_data, sim_data)

	## Setup constants
	def _buildConstants(self, raw_data, sim_data):
		mass_parameters = raw_data.mass_parameters
		self.__dict__.update(mass_parameters)

		self.cellDryMassFraction = 1. - self.cellWaterMassFraction
		self.avgCellToInitialCellConvFactor = np.exp(np.log(2) * self.avgCellCellCycleProgress)

	def _buildDependentConstants(self):
		self.avgCellDryMassInit = self.avgCellDryMass / self.avgCellToInitialCellConvFactor
		avgCellWaterMass = (self.avgCellDryMass / self.cellDryMassFraction) * self.cellWaterMassFraction
		self.avgCellWaterMassInit = avgCellWaterMass / self.avgCellToInitialCellConvFactor

	# Setup sub-masses
	def _buildSubMasses(self, raw_data, sim_data):
		self._doubling_time_vector = units.min * np.array([float(x['doublingTime'].asNumber(units.min)) for x in raw_data.dryMassComposition])

		# TODO: Use helper functions written for growthRateDependent parameters to make this better!
		dryMass = np.array([float(x['averageDryMass'].asNumber(units.fg)) for x in raw_data.dryMassComposition]) / self.avgCellToInitialCellConvFactor
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
	def _setAvgCellDryMass(self):
		doubling_time = self._clipTau_d(self._doubling_time)
		avgCellDryMass = units.fg * float(interpolate.splev(doubling_time.asNumber(units.min), self._dryMassParams))
		return avgCellDryMass

	# Set mass fractions based on growth rate
	def _setMassFraction(self):
		if type(self._doubling_time) != unum.Unum:
			raise Exception("Doubling time was not set!")

		D = {}
		dnaMassFraction = self._calculateGrowthRateDependentDnaMass(self._doubling_time) / self.avgCellDryMass
		dnaMassFraction.normalize()
		dnaMassFraction.checkNoUnit()
		D["dna"] = dnaMassFraction.asNumber()

		doubling_time = self._clipTau_d(self._doubling_time)
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

	def _setSubMass(self):
		D = {}
		for key, value in self.massFraction.iteritems():
			D[key + "Mass"] = value * self.avgCellDryMass

		D["rRna23SMass"] = D['rnaMass'] * self._rrna23s_mass_sub_fraction
		D["rRna16SMass"] = D['rnaMass'] * self._rrna16s_mass_sub_fraction
		D["rRna5SMass"] = D['rnaMass'] * self._rrna5s_mass_sub_fraction
		D["tRnaMass"] = D['rnaMass'] * self._trna_mass_sub_fraction
		D["mRnaMass"] = D['rnaMass'] * self._mrna_mass_sub_fraction

		return D

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
		growth_rate_unit = units.getUnit(raw_data.trna_growth_rates[0]['growth rate'])

		self._trna_growth_rates = growth_rate_unit * np.array([x['growth rate'].asNumber() for x in raw_data.trna_growth_rates])

		trna_ratio_to_16SrRNA_by_growth_rate = []
		for gr in self._trna_growth_rates: # This is a little crazy...
			trna_ratio_to_16SrRNA_by_growth_rate.append([x['ratio to 16SrRNA'] for x in getattr(raw_data, "trna_ratio_to_16SrRNA_" + str(gr.asNumber()).replace('.','p'))])
		self._trna_ratio_to_16SrRNA_by_growth_rate = np.array(trna_ratio_to_16SrRNA_by_growth_rate)

		self._trna_ids = [x['rna id'] for x in raw_data.trna_ratio_to_16SrRNA_0p4]

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

	def __init__(self, raw_data, sim_data):
		self._doubling_time = sim_data.doubling_time
		_loadTableIntoObjectGivenDoublingTime(self, raw_data.growthRateDependentParameters)
		# thingsToSet = [x for x in dir(self) if x[0] != '_']
		# for x in thingsToSet:
		# 	setattr(sim_data.constants, x, getattr(self, x))

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
