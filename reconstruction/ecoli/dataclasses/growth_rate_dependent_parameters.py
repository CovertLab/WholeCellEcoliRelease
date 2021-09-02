"""
SimulationData mass data
"""

from typing import Tuple

import numpy as np
from scipy import interpolate, stats
import unum

from wholecell.utils import fitting, units
import six


NORMAL_CRITICAL_MASS = 975 * units.fg
SLOW_GROWTH_FACTOR = 1.2  # adjustment for smaller cells


class Mass(object):
	""" Mass """

	def __init__(self, raw_data, sim_data):
		self._doubling_time = sim_data.doubling_time

		self._build_constants(raw_data, sim_data)
		self._build_submasses(raw_data, sim_data)
		self._build_CD_periods(raw_data, sim_data)

		self.avg_cell_dry_mass = self.get_avg_cell_dry_mass(self._doubling_time)
		self.mass_fractions = self.get_mass_fractions(self._doubling_time)
		self.avg_cell_component_masses = self.get_component_masses(self._doubling_time)

		self._build_dependent_constants()
		self._build_trna_data(raw_data, sim_data)

	## Setup constants
	def _build_constants(self, raw_data, sim_data):
		mass_parameters = raw_data.mass_parameters
		self.__dict__.update(mass_parameters)

		self.cell_dry_mass_fraction = 1. - self.cell_water_mass_fraction
		self.avg_cell_to_initial_cell_conversion_factor = np.exp(np.log(2) * self.avg_cell_cell_cycle_progress)

		self._cellDensity = sim_data.constants.cell_density

		self._glycogenFractions = raw_data.mass_fractions.glycogen_fractions
		self._mureinFractions = raw_data.mass_fractions.murein_fractions
		self._LPSFractions = raw_data.mass_fractions.LPS_fractions
		self._lipidFractions = raw_data.mass_fractions.lipid_fractions
		self._ionFractions = raw_data.mass_fractions.ion_fractions
		self._solubleFractions = raw_data.mass_fractions.soluble_fractions

		metIds = (
			[x["metaboliteId"] for x in self._glycogenFractions] +
			[x["metaboliteId"] for x in self._mureinFractions] +
			[x["metaboliteId"] for x in self._LPSFractions] +
			[x["metaboliteId"] for x in self._lipidFractions] +
			[x["metaboliteId"] for x in self._ionFractions] +
			[x["metaboliteId"] for x in self._solubleFractions] +
			["WATER[c]"]
			)
		mws = sim_data.getter.get_masses(metIds)
		self._mws = dict(zip(metIds, mws))

		self._metTargetIds = [x["Metabolite"] + "[c]" for x in raw_data.metabolite_concentrations]

	def _build_dependent_constants(self):
		self.avg_cell_dry_mass_init = self.avg_cell_dry_mass / self.avg_cell_to_initial_cell_conversion_factor
		avgCellWaterMass = (self.avg_cell_dry_mass / self.cell_dry_mass_fraction) * self.cell_water_mass_fraction
		self.avg_cell_water_mass_init = avgCellWaterMass / self.avg_cell_to_initial_cell_conversion_factor

	# Setup sub-masses
	def _build_submasses(self, raw_data, sim_data):
		self._doubling_time_vector = units.min * np.array([float(x['doublingTime'].asNumber(units.min)) for x in raw_data.dry_mass_composition])

		dryMass = np.array([
			float(x['averageDryMass'].asNumber(units.fg))
			for x in raw_data.dry_mass_composition
			])
		self._dryMassParams = linear_regression(
			self._doubling_time_vector.asNumber(units.min), 1. / dryMass)

		self._proteinMassFractionParams = self._getFitParameters(raw_data.dry_mass_composition, 'proteinMassFraction')
		self._rnaMassFractionParams = self._getFitParameters(raw_data.dry_mass_composition, 'rnaMassFraction')
		self._lipidMassFractionParams = self._getFitParameters(raw_data.dry_mass_composition, 'lipidMassFraction')
		self._lpsMassFractionParams = self._getFitParameters(raw_data.dry_mass_composition, 'lpsMassFraction')
		self._mureinMassFractionParams = self._getFitParameters(raw_data.dry_mass_composition, 'mureinMassFraction')
		self._glycogenMassFractionParams = self._getFitParameters(raw_data.dry_mass_composition, 'glycogenMassFraction')
		self._solublePoolMassFractionParams = self._getFitParameters(raw_data.dry_mass_composition, 'solublePoolMassFraction')
		self._inorganicIonMassFractionParams = self._getFitParameters(raw_data.dry_mass_composition, 'inorganicIonMassFraction')

		self.chromosome_sequence_mass = (
			sim_data.getter.get_mass(sim_data.molecule_ids.full_chromosome)
				/ sim_data.constants.n_avogadro
			).asUnit(units.g)

	def _getFitParameters(self, dryMassComposition, massFractionName):
		massFraction = np.array([float(x[massFractionName]) for x in dryMassComposition])
		x = self._doubling_time_vector.asNumber(units.min)[::-1]
		y = massFraction[::-1]
		massParams = interpolate.splrep(x, y)
		if np.sum(np.absolute(interpolate.splev(self._doubling_time_vector.asNumber(units.min), massParams) - massFraction)) / massFraction.size > 1.:
			raise Exception("Fitting {} with double exponential, residuals are huge!".format(massFractionName))
		return massParams

	def _build_CD_periods(self, raw_data, sim_data):
		self._c_period = sim_data.constants.c_period
		self._d_period = sim_data.constants.d_period

	# Set based on growth rate avgCellDryMass
	def get_avg_cell_dry_mass(self, doubling_time):
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

		mass = self.get_avg_cell_dry_mass(doubling_time) / self.cell_dry_mass_fraction
		critical_mass = min(mass * SLOW_GROWTH_FACTOR, NORMAL_CRITICAL_MASS)
		return critical_mass

	# Set mass fractions based on growth rate
	def get_mass_fractions(self, doubling_time):
		if type(doubling_time) != unum.Unum:
			raise Exception("Doubling time was not set!")

		D = {}
		dnaMassFraction = self._calculateGrowthRateDependentDnaMass(doubling_time) / self.get_avg_cell_dry_mass(doubling_time)
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


	def get_component_masses(self, doubling_time):
		D = {}
		massFraction = self.get_mass_fractions(doubling_time)
		for key, value in six.viewitems(massFraction):
			D[key + "Mass"] = value * self.get_avg_cell_dry_mass(doubling_time)

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

		avgCellDryMassInit = self.get_avg_cell_dry_mass(doubling_time) / self.avg_cell_to_initial_cell_conversion_factor
		avgCellWaterMassInit = (avgCellDryMassInit / self.cell_dry_mass_fraction) * self.cell_water_mass_fraction

		initWaterMass = avgCellWaterMassInit.asNumber(units.g)
		initDryMass = avgCellDryMassInit.asNumber(units.g)

		initCellMass = initWaterMass + initDryMass

		initCellVolume = initCellMass / self._cellDensity.asNumber(units.g / units.L) # L

		massFraction = self.get_mass_fractions(doubling_time)

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
		c_plus_d_period = self._c_period + self._d_period

		if doubling_time < self._d_period:
			raise Exception(
				"Can't have doubling time shorter than cytokinesis time!")

		doubling_time_unit = units.getUnit(doubling_time)

		# TODO: If you really care, this should be a loop.
		# It is optimized to run quickly over the range of T_d
		# and C and D periods that we have.
		return self.chromosome_sequence_mass * (1 +
												1 * (np.maximum(0. * doubling_time_unit, c_plus_d_period - doubling_time) / self._c_period) +
												2 * (np.maximum(0. * doubling_time_unit, c_plus_d_period - 2 * doubling_time) / self._c_period) +
												4 * (np.maximum(0. * doubling_time_unit, c_plus_d_period - 4 * doubling_time) / self._c_period)
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

	def _build_trna_data(self, raw_data, sim_data):
		growth_rate_unit = units.getUnit(raw_data.trna_data.trna_growth_rates[0]['growth rate'])

		self._trna_growth_rates = growth_rate_unit * np.array([x['growth rate'].asNumber() for x in raw_data.trna_data.trna_growth_rates])

		trna_ratio_to_16SrRNA_by_growth_rate = []
		for gr in self._trna_growth_rates: # This is a little crazy...
			trna_ratio_to_16SrRNA_by_growth_rate.append([x['ratio to 16SrRNA'] for x in getattr(raw_data.trna_data, "trna_ratio_to_16SrRNA_" + str(gr.asNumber()).replace('.','p'))])
		self._trna_ratio_to_16SrRNA_by_growth_rate = np.array(trna_ratio_to_16SrRNA_by_growth_rate)

		self._trna_ids = [x['rna id'] for x in raw_data.trna_data.trna_ratio_to_16SrRNA_0p4]

	def get_trna_distribution(self, doubling_time):
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
		_loadTableIntoObjectGivenDoublingTime(self, raw_data.growth_rate_dependent_parameters)
		self.ribosome_elongation_rate_params = _get_fit_parameters(raw_data.growth_rate_dependent_parameters, "ribosomeElongationRate")
		self.RNAP_elongation_rate_params = _get_fit_parameters(raw_data.growth_rate_dependent_parameters, "rnaPolymeraseElongationRate")
		self.RNAP_active_fraction_params = _get_fit_parameters(raw_data.growth_rate_dependent_parameters, "fractionActiveRnap")
		self.ribosome_active_fraction_params = _get_fit_parameters(raw_data.growth_rate_dependent_parameters, "fractionActiveRibosome")

		# ppGpp concentration by linear fit
		self.per_dry_mass_to_per_volume = sim_data.constants.cell_density * (1. - raw_data.mass_parameters['cell_water_mass_fraction'])
		doubling_time = _loadRow('doublingTime', raw_data.growth_rate_dependent_parameters)
		ppgpp_conc = _loadRow('ppGpp_conc', raw_data.growth_rate_dependent_parameters) * self.per_dry_mass_to_per_volume
		self._ppGpp_concentration = _get_linearized_fit(doubling_time, ppgpp_conc)

	def get_ribosome_elongation_rate(self, doubling_time):
		return _useFitParameters(doubling_time, **self.ribosome_elongation_rate_params)

	def get_rnap_elongation_rate(self, doubling_time):
		return _useFitParameters(doubling_time, **self.RNAP_elongation_rate_params)

	def get_fraction_active_rnap(self, doubling_time):
		return _useFitParameters(doubling_time, **self.RNAP_active_fraction_params)

	def get_fraction_active_ribosome(self, doubling_time):
		return _useFitParameters(doubling_time, **self.ribosome_active_fraction_params)

	def get_ppGpp_conc(self, doubling_time):
		return _use_linearized_fit(doubling_time, self._ppGpp_concentration)

def _get_fit_parameters(list_of_dicts, key):
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

def _get_linearized_fit(x, y, **kwargs):
	if units.hasUnit(x):
		x_units = units.getUnit(x)
		x = x.asNumber(x_units)
	else:
		x_units = None
	if units.hasUnit(y):
		y_units = units.getUnit(y)
		y = y.asNumber(y_units)
	else:
		y_units = 1.

	return x_units, y_units, fitting.fit_linearized_transforms(x, y, **kwargs)

def _use_linearized_fit(x, params):
	x_units, y_units, fit_params = params
	x_unitless = x if x_units is None else x.asNumber(x_units)
	return y_units * fitting.interpolate_linearized_fit(x_unitless, *fit_params)

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
	if dtype == int:
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
		fitParameters = _get_fit_parameters(list_of_dicts, key)
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
