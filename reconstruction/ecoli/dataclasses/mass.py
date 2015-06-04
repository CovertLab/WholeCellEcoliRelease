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
		self.massFraction = self._setMassFraction()
		self.subMass = self._setSubMass()

		self._buildDependentConstants()

		self._buildTrnaData(raw_data, sim_data)


	## Setup constants
	def _buildConstants(self, raw_data, sim_data):
		mass_parameters = raw_data.mass_parameters
		self.__dict__.update(mass_parameters)

		self.cellDryMassFraction = 1. - self.cellWaterMassFraction
		self.avgCellToInitalCellConvFactor = np.exp(np.log(2) * self.avgCellCellCycleProgress)

	def _buildDependentConstants(self):
		self.avgCellDryMassInit = self.avgCellDryMass / self.avgCellToInitalCellConvFactor
		avgCellWaterMass = (self.avgCellDryMass / self.cellDryMassFraction) * self.cellWaterMassFraction
		self.avgCellWaterMassInit = avgCellWaterMass / self.avgCellToInitalCellConvFactor

	# Setup sub-masses
	def _buildSubMasses(self, raw_data, sim_data):
		self._doubling_time_vector = units.min * np.array([float(x['doublingTime'].asNumber(units.min)) for x in raw_data.dryMassComposition])

		dryMass = np.array([float(x['averageDryMass'].asNumber(units.fg)) for x in raw_data.dryMassComposition]) / self.avgCellToInitalCellConvFactor
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
		self.c_period = sim_data.constants.c_period
		self.d_period = sim_data.constants.d_period

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
		assert np.absolute(np.sum([x for x in D.itervalues()]) - 1.) < 1e3
		return D

	def _setSubMass(self):
		D = {}
		for key, value in self.massFraction.iteritems():
			D[key + "Mass"] = value * self.avgCellDryMass

		D["rRna23SMass"] = D['rnaMass'] * self.rrna23s_mass_sub_fraction
		D["rRna16SMass"] = D['rnaMass'] * self.rrna16s_mass_sub_fraction
		D["rRna5SMass"] = D['rnaMass'] * self.rrna5s_mass_sub_fraction
		D["tRnaMass"] = D['rnaMass'] * self.trna_mass_sub_fraction
		D["mRnaMass"] = D['rnaMass'] * self.mrna_mass_sub_fraction

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
