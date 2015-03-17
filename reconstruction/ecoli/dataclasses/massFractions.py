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

class MassFractions(object):
	""" MassFractions """

	def __init__(self, raw_data, sim_data):
		self._buildMassFractions(raw_data, sim_data)

	def _buildMassFractions(self, raw_data, sim_data):
		self.tau_d = units.min * np.array([float(x['doublingTime'].asNumber(units.min)) for x in raw_data.dryMassComposition])

		avgToBeginningConvFactor = raw_data.parameters['avgCellToInitalCellConvFactor']
		self._dryMass = np.array([float(x['averageDryMass'].asNumber(units.fg)) for x in raw_data.dryMassComposition]) / avgToBeginningConvFactor
		self._proteinMass = self._dryMass * np.array([float(x['proteinMassFraction']) for x in raw_data.dryMassComposition])
		self._rnaMass = self._dryMass * np.array([float(x['rnaMassFraction']) for x in raw_data.dryMassComposition])
		self._dnaMass = self._dryMass * np.array([float(x['dnaMassFraction']) for x in raw_data.dryMassComposition])

		# We are assuming these are constant over all growth rates
		# (probably not be true...)
		self.RRNA23S_MASS_SUB_FRACTION = raw_data.parameters['rrna23s_mass_sub_fraction']
		self.RRNA16S_MASS_SUB_FRACTION = raw_data.parameters['rrna16s_mass_sub_fraction']
		self.RRNA5S_MASS_SUB_FRACTION = raw_data.parameters['rrna5s_mass_sub_fraction']
		self.TRNA_MASS_SUB_FRACTION = raw_data.parameters['trna_mass_sub_fraction']
		self.MRNA_MASS_SUB_FRACTION = raw_data.parameters['mrna_mass_sub_fraction']

		self.dryMassParams, _ = curve_fit(self._exp2, self.tau_d.asNumber(units.min), self._dryMass, p0 = (0, 0, 0, 0))
		self.proteinMassParams, _ = curve_fit(self._exp2, self.tau_d.asNumber(units.min), self._proteinMass, p0 = (0, 0, 0, 0))
		self.rnaMassParams, _ = curve_fit(self._exp2, self.tau_d.asNumber(units.min), self._rnaMass, p0 = (0, 0, 0, 0))
		self.dnaMassParams, _ = curve_fit(self._exp2, self.tau_d.asNumber(units.min), self._dnaMass, p0 = (0, 0, 0, 0))

		self.chromMass = self._chromMass(raw_data, sim_data)

		self.C_PERIOD = raw_data.parameters['c_period']
		self.D_PERIOD = raw_data.parameters['d_period']
		self.CD_PERIOD = self.C_PERIOD + self.D_PERIOD

	def massFractions(self, tau_d):
		"""
		Given an input doubling time in minutes, output mass fractions in fg
		"""

		D = {}
		D["dnaMass"] = self._calculateDnaMass(tau_d)

		tau_d = self._clipTau_d(tau_d)

		D["proteinMass"] = units.fg * self._exp2(tau_d.asNumber(units.min), *self.proteinMassParams)
		D["rnaMass"] = units.fg * self._exp2(tau_d.asNumber(units.min), *self.rnaMassParams)
		D["rRna23SMass"] = D["rnaMass"] * self.RRNA23S_MASS_SUB_FRACTION
		D["rRna16SMass"] = D["rnaMass"] * self.RRNA16S_MASS_SUB_FRACTION
		D["rRna5SMass"] = D["rnaMass"] * self.RRNA5S_MASS_SUB_FRACTION
		D["tRnaMass"] = D["rnaMass"] * self.TRNA_MASS_SUB_FRACTION
		D["mRnaMass"] = D["rnaMass"] * self.MRNA_MASS_SUB_FRACTION

		return D

	def _calculateDnaMass(self, tau_d):
		if tau_d < self.D_PERIOD:
			raise Exception, "Can't have doubling time shorter than cytokinesis time!"

		tau_d_unit = units.getUnit(tau_d)

		# TODO: If you really care, this should be a loop.
		# It is optimized to run quickly over the range of T_d
		# and C and D periods that we have.
		return self.chromMass * (1 +
			1 * (np.maximum(0. * tau_d_unit, self.CD_PERIOD - tau_d) / self.C_PERIOD) +
			2 * (np.maximum(0. * tau_d_unit, self.CD_PERIOD - 2 * tau_d) / self.C_PERIOD) +
			4 * (np.maximum(0. * tau_d_unit, self.CD_PERIOD - 4 * tau_d) / self.C_PERIOD)
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
			tau_d[tau_d > max(self.tau_d)] = max(self.tau_d)
			tau_d[tau_d < min(self.tau_d)] = min(self.tau_d)
		else:
			if tau_d > max(self.tau_d):
				tau_d = max(self.tau_d)
			elif tau_d < min(self.tau_d):
				tau_d = min(self.tau_d)
		return tau_d

	def _exp2(self, x, a, b, c, d):
		return a * np.exp(b * x) + c * np.exp(d * x)