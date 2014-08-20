#!/usr/bin/env python

import numpy as np
from scipy.optimize import curve_fit

def exp2(x, a, b, c, d):
	return a * np.exp(b * x) + c * np.exp(d * x)

class GrowthData(object):

	def __init__(self):
		self.tau_d = np.array([100., 60., 40., 30., 24.])
		self.mu = np.array([0.6, 1, 1.5, 2, 2.5])

		avgToBeginningConvFactor = 1.36;
		self._dryMass = np.array([148., 258., 433., 641., 865.]) / avgToBeginningConvFactor
		self._proteinMass = np.array([100., 156., 234., 340., 450.]) / avgToBeginningConvFactor
		self._rnaMass = np.array([20., 39., 77., 132., 211.]) / avgToBeginningConvFactor
		self._dnaMass = np.array([7.6, 9, 11.3, 14.4, 18.3]) / avgToBeginningConvFactor

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

	def massFractions(self, tau_d):
		"""
		Given an input doubling time in minutes, output mass fractions in fg
		"""

		# Clip values to be in the range that we have data for
		if hasattr(tau_d, "dtype"):
			tau_d[tau_d > self.tau_d.max()] = self.tau_d.max()
			tau_d[tau_d < self.tau_d.min()] = self.tau_d.min()
		else:
			if tau_d > self.tau_d.max():
				tau_d = self.tau_d.max()
			elif tau_d < self.tau_d.min():
				tau_d = self.tau_d.min()

		D = {}
		D["dryMass"] = exp2(tau_d, *self.dryMassParams)
		D["proteinMass"] = exp2(tau_d, *self.proteinMassParams)
		D["rnaMass"] = exp2(tau_d, *self.rnaMassParams)
		D["dnaMass"] = exp2(tau_d, *self.dnaMassParams)
		D["rRna23SMass"] = D["rnaMass"] * self.RRNA23S_MASS_SUB_FRACTION
		D["rRna16SMass"] = D["rnaMass"] * self.RRNA16S_MASS_SUB_FRACTION
		D["rRna5SMass"] = D["rnaMass"] * self.RRNA5S_MASS_SUB_FRACTION
		D["tRnaMass"] = D["rnaMass"] * self.TRNA_MASS_SUB_FRACTION
		D["mRnaMass"] = D["rnaMass"] * self.MRNA_MASS_SUB_FRACTION

		return D