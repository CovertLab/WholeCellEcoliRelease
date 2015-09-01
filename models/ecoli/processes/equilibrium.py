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
import theano.tensor as T
import theano

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

		self.nAvogadro = kb.constants.nAvogadro.asNumber(1 / units.mol)
		self.cellDensity = kb.constants.cellDensity.asNumber(units.g / units.L)

		# Create matrices and vectors

		self.stoichMatrix = kb.process.equilibrium.stoichMatrix().astype(np.int64)
		self.ratesFwd = kb.process.equilibrium.ratesFwd
		self.ratesRev = kb.process.equilibrium.ratesRev

		self._makeDerivative()
		self._makeMatrices()

		# Build views

		moleculeNames = kb.process.equilibrium.moleculeNames

		self.molecules = self.bulkMoleculesView(moleculeNames)


	def calculateRequest(self):
		moleculeCounts = self.molecules.total()

		cellMass = (self.readFromListener("Mass", "cellMass") * units.fg).asNumber(units.g)
		cellVolume = cellMass / self.cellDensity

		y_init = moleculeCounts / (cellVolume * self.nAvogadro)
		y = scipy.integrate.odeint(self._derivatives, y_init, t = [0, 1e4], Dfun = self._derivativesJacobian)

		if np.any(y[-1, :] * (cellVolume * self.nAvogadro) <= -1):
			raise Exception, "Have negative values -- probably due to numerical instability"
		if np.linalg.norm(self._derivatives(y[-1, :], 0), np.inf) * (cellVolume * self.nAvogadro) > 1:
			raise Exception, "Didn't reach steady state"
		y[y < 0] = 0
		yMolecules = y * (cellVolume * self.nAvogadro)

		dYMolecules = yMolecules[-1, :] - yMolecules[0, :]
		rxnFluxes = np.round(np.dot(self.metsToRxnFluxes, dYMolecules))
		rxnFluxesN = -1. * (rxnFluxes < 0) * rxnFluxes
		rxnFluxesP =  1. * (rxnFluxes > 0) * rxnFluxes

		self.req = np.dot(self.Rp, rxnFluxesP) + np.dot(self.Pp, rxnFluxesN)
		self.rxnFluxes = rxnFluxes

		self.molecules.requestIs(self.req)

	def evolveState(self):
		moleculeCounts = self.molecules.counts()

		rxnFluxes = self.rxnFluxes.copy()

		# Kill "bad" reactions where we have insufficient metabolites
		for badMetIdx in np.where(self.req > moleculeCounts)[0]:
			rxnFluxes[self.stoichMatrix[badMetIdx, :]!= 0] = 0

		assert(np.all(moleculeCounts + np.dot(self.stoichMatrix, rxnFluxes) >= 0))

		self.molecules.countsInc(
			np.dot(self.stoichMatrix, rxnFluxes)
			)

	def _makeDerivative(self):

		S = self.stoichMatrix

		y = T.dvector()
		dy = [0 * y[0] for _ in xrange(S.shape[0])]
		for colIdx in xrange(S.shape[1]):
			negIdxs = np.where(S[:, colIdx] < 0)[0]
			posIdxs = np.where(S[:, colIdx] > 0)[0]

			reactantFlux = self.ratesFwd[colIdx]
			for negIdx in negIdxs:
				reactantFlux *= (y[negIdx] ** (-1 * S[negIdx, colIdx]))

			productFlux = self.ratesRev[colIdx]
			for posIdx in posIdxs:
				productFlux *=  (y[posIdx] ** ( 1 * S[posIdx, colIdx]))

			fluxForNegIdxs = (-1. * reactantFlux) + (1. * productFlux)
			fluxForPosIdxs = ( 1. * reactantFlux) - (1. * productFlux)

			for thisIdx in negIdxs:
				dy[thisIdx] += fluxForNegIdxs
			for thisIdx in posIdxs:
				dy[thisIdx] += fluxForPosIdxs

		t = T.dscalar()

		J = [T.grad(dy[i], y) for i in xrange(len(dy))]

		self._derivativesJacobian = theano.function([y, t], T.stack(*J), on_unused_input = "ignore")
		self._derivatives = theano.function([y, t], T.stack(*dy), on_unused_input = "ignore")

	def _makeMatrices(self):
		EPS = 1e-9

		S = self.stoichMatrix
		S1 = np.zeros_like(S)
		S1[S < -1 * EPS] = -1
		S1[S > EPS] = 1

		Rp =  1. * (S1 < 0)
		Pp =  1. * (S1 > 0)

		self.Rp = Rp
		self.Pp = Pp

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
