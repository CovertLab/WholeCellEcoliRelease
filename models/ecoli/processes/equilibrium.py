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
		self.Rp = kb.process.equilibrium.Rp
		self.Pp = kb.process.equilibrium.Pp
		self.derivatives = kb.process.equilibrium.derivatives
		self.derivativesJacobian = kb.process.equilibrium.derivativesJacobian
		self.metsToRxnFluxes = kb.process.equilibrium.metsToRxnFluxes

		self.fluxesAndMoleculesToSS = kb.process.equilibrium.fluxesAndMoleculesToSS

		# Build views

		moleculeNames = kb.process.equilibrium.moleculeNames

		self.molecules = self.bulkMoleculesView(moleculeNames)


	def calculateRequest(self):
		moleculeCounts = self.molecules.total()

		cellMass = (self.readFromListener("Mass", "cellMass") * units.fg).asNumber(units.g)
		cellVolume = cellMass / self.cellDensity

		self.rxnFluxes, self.req = self.fluxesAndMoleculesToSS(moleculeCounts, cellVolume, self.nAvogadro)

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



