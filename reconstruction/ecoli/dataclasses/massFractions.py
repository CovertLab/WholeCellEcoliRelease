"""
SimulationData mass fraction data

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 03/13/2015
"""

from __future__ import division

import numpy as np
from scipy.optimize import curve_fit
from wholecell.utils import units
import unum

class MassFractions(object):
	""" MassFractions """

	def __init__(self, raw_data, sim_data):
		self._buildMassFractions(raw_data, sim_data)
		self._buildTrnaData(raw_data, sim_data)

	def _buildMassFractions(self, raw_data, sim_data):
		self._tau_d = units.min * np.array([float(x['doublingTime'].asNumber(units.min)) for x in raw_data.dryMassComposition])

		avgToBeginningConvFactor = raw_data.parameters['avgCellToInitalCellConvFactor']
		self._dryMass = np.array([float(x['averageDryMass'].asNumber(units.fg)) for x in raw_data.dryMassComposition]) / avgToBeginningConvFactor
		self._proteinMass = self._dryMass * np.array([float(x['proteinMassFraction']) for x in raw_data.dryMassComposition])
		self._rnaMass = self._dryMass * np.array([float(x['rnaMassFraction']) for x in raw_data.dryMassComposition])
		self._dnaMass = self._dryMass * np.array([float(x['dnaMassFraction']) for x in raw_data.dryMassComposition])

		# We are assuming these are constant over all growth rates
		# (probably not be true...)
		self._RRNA23S_MASS_SUB_FRACTION = raw_data.parameters['rrna23s_mass_sub_fraction']
		self._RRNA16S_MASS_SUB_FRACTION = raw_data.parameters['rrna16s_mass_sub_fraction']
		self._RRNA5S_MASS_SUB_FRACTION = raw_data.parameters['rrna5s_mass_sub_fraction']
		self._TRNA_MASS_SUB_FRACTION = raw_data.parameters['trna_mass_sub_fraction']
		self._MRNA_MASS_SUB_FRACTION = raw_data.parameters['mrna_mass_sub_fraction']

		self._dryMassParams, _ = curve_fit(self._exp2, self._tau_d.asNumber(units.min), self._dryMass, p0 = (0, 0, 0, 0))
		self._proteinMassParams, _ = curve_fit(self._exp2, self._tau_d.asNumber(units.min), self._proteinMass, p0 = (0, 0, 0, 0))
		self._rnaMassParams, _ = curve_fit(self._exp2, self._tau_d.asNumber(units.min), self._rnaMass, p0 = (0, 0, 0, 0))
		self._dnaMassParams, _ = curve_fit(self._exp2, self._tau_d.asNumber(units.min), self._dnaMass, p0 = (0, 0, 0, 0))

		self.chromMass = self._chromMass(raw_data, sim_data)

		self._C_PERIOD = raw_data.parameters['c_period']
		self._D_PERIOD = raw_data.parameters['d_period']
		self._CD_PERIOD = self._C_PERIOD + self._D_PERIOD

	def massFractions(self, tau_d):
		"""
		Given an input doubling time in minutes, output mass fractions in fg
		"""

		D = {}
		D["dnaMass"] = self._calculateDnaMass(tau_d)

		tau_d = self._clipTau_d(tau_d)

		D["proteinMass"] = units.fg * self._exp2(tau_d.asNumber(units.min), *self._proteinMassParams)
		D["rnaMass"] = units.fg * self._exp2(tau_d.asNumber(units.min), *self._rnaMassParams)
		D["rRna23SMass"] = D["rnaMass"] * self._RRNA23S_MASS_SUB_FRACTION
		D["rRna16SMass"] = D["rnaMass"] * self._RRNA16S_MASS_SUB_FRACTION
		D["rRna5SMass"] = D["rnaMass"] * self._RRNA5S_MASS_SUB_FRACTION
		D["tRnaMass"] = D["rnaMass"] * self._TRNA_MASS_SUB_FRACTION
		D["mRnaMass"] = D["rnaMass"] * self._MRNA_MASS_SUB_FRACTION

		return D

	def _calculateDnaMass(self, tau_d):
		if tau_d < self._D_PERIOD:
			raise Exception, "Can't have doubling time shorter than cytokinesis time!"

		tau_d_unit = units.getUnit(tau_d)

		# TODO: If you really care, this should be a loop.
		# It is optimized to run quickly over the range of T_d
		# and C and D periods that we have.
		return self.chromMass * (1 +
			1 * (np.maximum(0. * tau_d_unit, self._CD_PERIOD - tau_d) / self._C_PERIOD) +
			2 * (np.maximum(0. * tau_d_unit, self._CD_PERIOD - 2 * tau_d) / self._C_PERIOD) +
			4 * (np.maximum(0. * tau_d_unit, self._CD_PERIOD - 4 * tau_d) / self._C_PERIOD)
			)

	def _chromMass(self, raw_data, sim_data):
		dntCounts = np.array([
			raw_data.genome_sequence.count('A') + raw_data.genome_sequence.count('T'),
			raw_data.genome_sequence.count('C') + raw_data.genome_sequence.count('G'),
			raw_data.genome_sequence.count('G') + raw_data.genome_sequence.count('C'),
			raw_data.genome_sequence.count('T') + raw_data.genome_sequence.count('A')
		])

		dntMasses = (sim_data.getter.getMass(sim_data.moleculeGroups.polymerizedDNT_IDs) / raw_data.constants['nAvogadro']).asUnit(units.g)
		chromMass = units.dot(dntCounts, dntMasses)
		return chromMass

	def _clipTau_d(self, tau_d):
		# Clip values to be in the range that we have data for
		if hasattr(tau_d, "dtype"):
			tau_d[tau_d > max(self._tau_d)] = max(self._tau_d)
			tau_d[tau_d < min(self._tau_d)] = min(self._tau_d)
		else:
			if tau_d > max(self._tau_d):
				tau_d = max(self._tau_d)
			elif tau_d < min(self._tau_d):
				tau_d = min(self._tau_d)
		return tau_d

	def _exp2(self, x, a, b, c, d):
		return a * np.exp(b * x) + c * np.exp(d * x)

	def _buildTrnaData(self, raw_data, sim_data):
		growth_rate_unit = units.getUnit(raw_data.trna_growth_rates[0]['growth rate'])

		self._trna_growth_rates = growth_rate_unit * np.array([x['growth rate'].asNumber() for x in raw_data.trna_growth_rates])

		trna_ratio_to_16SrRNA_by_growth_rate = []
		for gr in self._trna_growth_rates: # This is a little crazy...
			trna_ratio_to_16SrRNA_by_growth_rate.append([x['ratio to 16SrRNA'] for x in getattr(raw_data, "trna_ratio_to_16SrRNA_" + str(gr.asNumber()).replace('.','p'))])
		self._trna_ratio_to_16SrRNA_by_growth_rate = np.array(trna_ratio_to_16SrRNA_by_growth_rate)

		self._trna_ids = [x['rna id'] for x in raw_data.trna_ratio_to_16SrRNA_0p4]

	def getTrnaAbundanceAtGrowthRate(self, tau_d):
		assert type(tau_d) == unum.Unum
		assert type(tau_d.asNumber()) == float
		growth_rate = 1 / tau_d
		growth_rate = growth_rate.asNumber(1/units.h)

		from scipy.interpolate import interp1d
		trna_abundance_interpolation_functions = [interp1d(self._trna_growth_rates.asNumber(1/units.h), self._trna_ratio_to_16SrRNA_by_growth_rate[:,i]) for i in range(self._trna_ratio_to_16SrRNA_by_growth_rate.shape[1])]

		abundance = np.zeros(len(self._trna_ids), dtype = [('id','a50'),('molar_ratio_to_16SrRNA', np.float64)])
		abundance['id'] = self._trna_ids
		abundance['molar_ratio_to_16SrRNA'] = [x(growth_rate) for x in trna_abundance_interpolation_functions]
		return abundance





