#!/usr/bin/env python

"""
Equilibrium

Equilibrium binding sub-model

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 8/21/2015
"""

from __future__ import division

import numpy as np
import scipy.integrate
from wholecell.utils import units

import wholecell.processes.process

class Equilibrium(wholecell.processes.process.Process):
	""" Equilibrium """

	_name = "Equilibrium"

	# Constructor
	def __init__(self):

		super(Equilibrium, self).__init__()


	# Construct object graph
	def initialize(self, sim, kb):
		super(Equilibrium, self).initialize(sim, kb)

		# Create matrices and vectors

		self.stoichMatrix = kb.process.equilibrium.stoichMatrix().astype(np.int64)
		self.ratesFwd = kb.process.equilibrium.ratesFwd
		self.ratesRev = kb.process.equilibrium.ratesRev

		self._makeMatrices()

		self.nAvogadro = kb.constants.nAvogadro.asNumber(1 / units.mol)
		self.cellDensity = kb.constants.cellDensity.asNumber(units.g / units.L)

		# Build views

		moleculeNames = kb.process.equilibrium.moleculeNames

		self.molecules = self.bulkMoleculesView(moleculeNames)


	def calculateRequest(self):
		moleculeCounts = self.molecules.total()

		cellMass = (self.readFromListener("Mass", "cellMass") * units.fg).asNumber(units.g)
		cellVolume = cellMass / self.cellDensity

		y = scipy.integrate.odeint(self._derivatives, moleculeCounts / (self.nAvogadro * cellVolume), t = [0, 1000])
		yMolecules = y * (self.nAvogadro * cellVolume)

		dYMolecules = yMolecules[-1, :] - yMolecules[0, :]
		rxnFluxes = np.round(np.dot(self.metsToRxnFluxes, dYMolecules))
		rxnFluxesN = -1. * (rxnFluxes < 0) * rxnFluxes
		rxnFluxesP =  1. * (rxnFluxes > 0) * rxnFluxes

		self.req = np.dot(self.Rp, rxnFluxesP) + np.dot(self.Pp, rxnFluxesN)
		self.rxnFluxes = rxnFluxes

		self.molecules.requestIs(self.req)

	def evolveState(self):
		moleculeCounts = self.molecules.counts()

		# Don't do anything this time step if you don't get all the molecules requested
		if not np.all(self.req.astype(np.int) == moleculeCounts):
			return

		assert(np.all(moleculeCounts + np.dot(self.stoichMatrix, self.rxnFluxes) >= 0))

		self.molecules.countsInc(
			np.dot(self.stoichMatrix, self.rxnFluxes)
			)

	def _makeMatrices(self):
		EPS = 1e-9

		S = self.stoichMatrix
		S1 = np.zeros_like(S)
		S1[S < -1 * EPS] = -1
		S1[S > EPS] = 1

		R = (-1 * (S < 0) * S).T
		P = (	  (S > 0) * S).T

		Rp =  1. * (S1 < 0)
		Rn = -1. * (S1 < 0)
		Pp =  1. * (S1 > 0)
		Pn = -1. * (S1 > 0)

		self.S1 = S1
		self.R = scipy.sparse.csr_matrix(R)
		self.P = scipy.sparse.csr_matrix(P)
		self.Rp = Rp
		self.Rn = Rn
		self.Pp = Pp
		self.Pn = Pn

		metsToRxnFluxes = self.stoichMatrix.copy()

		metsToRxnFluxes[(np.abs(metsToRxnFluxes) > EPS).sum(axis = 1) > 1, : ] = 0
		for colIdx in xrange(metsToRxnFluxes.shape[1]):
			try:
				firstNonZeroIdx = np.where(np.abs(metsToRxnFluxes[:, colIdx]) > EPS)[0][0]
			except IndexError:
				raise Exception, "Column %d of S matrix not linearly independent!" % colIdx
			metsToRxnFluxes[:firstNonZeroIdx, colIdx] = 0
			metsToRxnFluxes[(firstNonZeroIdx + 1):, colIdx] = 0

		self.metsToRxnFluxes = metsToRxnFluxes.T

	def _derivatives(self, v, t):
		v[v < 0] = 0
		vLog = np.log(v)
		fluxReactant = self.ratesFwd * np.exp(self.R.dot(vLog))
		fluxProduct = self.ratesRev * np.exp(self.P.dot(vLog))

		return (
		 np.dot(self.Rn, fluxReactant) + np.dot(self.Rp, fluxProduct) +
		 np.dot(self.Pp, fluxReactant) + np.dot(self.Pn, fluxProduct)
		)