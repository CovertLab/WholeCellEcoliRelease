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
		self._makeMatrices()

		# Build views

		moleculeNames = kb.process.equilibrium.moleculeNames

		self.molecules = self.bulkMoleculesView(moleculeNames)


	def calculateRequest(self):
		moleculeCounts = self.molecules.total()

		self.molecules.requestIs(np.zeros_like(moleculeCounts))


	def evolveState(self):
		moleculeCounts = self.molecules.counts()

		self.molecules.countsIs(moleculeCounts)

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
		self.R = R
		self.P = P
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

		self.metsToRxnFluxes = metsToRxnFluxes

	def _derivatives(self, v, t):

		vLog = np.log(v)
		fluxReactant = self.ratesFwd * np.exp(np.dot(self.R, vLog))
		fluxProduct = self.ratesRev * np.exp(np.dot(self.P, vLog))
		return (
		 np.dot(self.Rn, fluxReactant) + np.dot(self.Rp, fluxProduct) +
		 np.dot(self.Pp, fluxReactant) + np.dot(self.Pn, fluxProduct)
		)