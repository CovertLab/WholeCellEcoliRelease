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
	def initialize(self, sim, sim_data):
		super(Equilibrium, self).initialize(sim, sim_data)

		self.nAvogadro = sim_data.constants.nAvogadro.asNumber(1 / units.mol)
		self.cellDensity = sim_data.constants.cellDensity.asNumber(units.g / units.L)

		# Create matrices and vectors

		self.stoichMatrix = sim_data.process.equilibrium.stoichMatrix().astype(np.int64)
		self.Rp = sim_data.process.equilibrium.Rp
		self.Pp = sim_data.process.equilibrium.Pp
		self.derivatives = sim_data.process.equilibrium.derivatives
		self.derivativesJacobian = sim_data.process.equilibrium.derivativesJacobian
		self.metsToRxnFluxes = sim_data.process.equilibrium.metsToRxnFluxes

		self.fluxesAndMoleculesToSS = sim_data.process.equilibrium.fluxesAndMoleculesToSS

		# Build views

		moleculeNames = sim_data.process.equilibrium.moleculeNames

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



