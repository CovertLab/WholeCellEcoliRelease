#!/usr/bin/env python

import numpy as np
from scipy.optimize import curve_fit
from wholecell.utils import units

def exp2(x, a, b, c, d):
	return a * np.exp(b * x) + c * np.exp(d * x)

class GrowthData(object):

	def __init__(self, kb):
		self.tau_d = np.array(kb.cellDryMassComposition["doublingTime"].asNumber(units.min))

		avgToBeginningConvFactor = kb.avgCellToInitalCellConvFactor
		self._dryMass = np.array([148., 258., 433., 641., 865.]) / avgToBeginningConvFactor # TOKB
		self._proteinMass = self._dryMass * kb.cellDryMassComposition["proteinMassFraction"]
		self._rnaMass = self._dryMass * kb.cellDryMassComposition["rnaMassFraction"]
		self._dnaMass = self._dryMass * kb.cellDryMassComposition["dnaMassFraction"]

		# We are assuming these are constant over all growth rates
		# (probably not be true...)
		self.RRNA23S_MASS_SUB_FRACTION = 0.525 # This is the fraction of RNA that is 23S rRNA
		self.RRNA16S_MASS_SUB_FRACTION = 0.271 # This is the fraction of RNA that is 16S rRNA
		self.RRNA5S_MASS_SUB_FRACTION = 0.017 # This is the fraction of RNA that is 5S rRNA
		self.TRNA_MASS_SUB_FRACTION = 0.146 # This is the fraction of RNA that is tRNA
		self.MRNA_MASS_SUB_FRACTION = 0.041 # This is the fraction of RNA that is mRNA

		self.dryMassParams, _ = curve_fit(exp2, self.tau_d, self._dryMass, p0 = (0, 0, 0, 0))
		self.proteinMassParams, _ = curve_fit(exp2, self.tau_d, self._proteinMass, p0 = (0, 0, 0, 0))
		self.rnaMassParams, _ = curve_fit(exp2, self.tau_d, self._rnaMass, p0 = (0, 0, 0, 0))
		self.dnaMassParams, _ = curve_fit(exp2, self.tau_d, self._dnaMass, p0 = (0, 0, 0, 0))

		self.chromMass = self._chromMass(kb)

		self.C_PERIOD = 40. # TOKB. [minutes]
		self.D_PERIOD = 20. # TOKB. [minutes]
		self.CD_PERIOD = self.C_PERIOD + self.D_PERIOD

	def _chromMass(self, kb):
		dntCounts = np.array([
			kb.genome_A_count + kb.genome_T_count,
			kb.genome_C_count + kb.genome_G_count,
			kb.genome_G_count + kb.genome_C_count,
			kb.genome_T_count + kb.genome_A_count
		])

		dntMasses = (kb.getMass(kb.moleculeGroups.polymerizedDNT_IDs) / kb.nAvogadro).asUnit(units.g)

		chromMass = units.dot(dntCounts, dntMasses)
		return chromMass

	def _clipTau_d(self, tau_d):
		# Clip values to be in the range that we have data for
		if hasattr(tau_d, "dtype"):
			tau_d[tau_d > self.tau_d.max()] = self.tau_d.max()
			tau_d[tau_d < self.tau_d.min()] = self.tau_d.min()
		else:
			if tau_d > self.tau_d.max():
				tau_d = self.tau_d.max()
			elif tau_d < self.tau_d.min():
				tau_d = self.tau_d.min()
		return tau_d


	def dnaMass(self, tau_d):
		if tau_d < self.D_PERIOD:
			raise Exception, "Can't have doubling time shorter than cytokinesis time!"

		# TODO: If you really care, this should be a loop.
		# It is optimized to run quickly over the range of T_d
		# and C and D periods that we have.
		return self.chromMass * (1 +
			1 * (np.maximum(0., self.CD_PERIOD - tau_d) / self.C_PERIOD) +
			2 * (np.maximum(0., self.CD_PERIOD - 2 * tau_d) / self.C_PERIOD) +
			4 * (np.maximum(0., self.CD_PERIOD - 4 * tau_d) / self.C_PERIOD)
			)


	def massFractions(self, tau_d):
		"""
		Given an input doubling time in minutes, output mass fractions in fg
		"""

		D = {}
		D["dnaMass"] = self.dnaMass(tau_d)

		tau_d = self._clipTau_d(tau_d)

		D["proteinMass"] = units.fg * exp2(tau_d, *self.proteinMassParams)
		D["rnaMass"] = units.fg * exp2(tau_d, *self.rnaMassParams)
		D["rRna23SMass"] = D["rnaMass"] * self.RRNA23S_MASS_SUB_FRACTION
		D["rRna16SMass"] = D["rnaMass"] * self.RRNA16S_MASS_SUB_FRACTION
		D["rRna5SMass"] = D["rnaMass"] * self.RRNA5S_MASS_SUB_FRACTION
		D["tRnaMass"] = D["rnaMass"] * self.TRNA_MASS_SUB_FRACTION
		D["mRnaMass"] = D["rnaMass"] * self.MRNA_MASS_SUB_FRACTION

		return D