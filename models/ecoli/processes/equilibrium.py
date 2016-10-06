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

		# If we didn't get allocated all the molecules we need, make do with what we have
		# (decrease reaction fluxes so that they make use of what we have, but not more)
		insufficientMetaboliteIdxs = np.where(self.req > moleculeCounts)[0]
		for insufficientMetaboliteIdx in insufficientMetaboliteIdxs:
			rxnPosIdxs = np.where(np.logical_and(self.stoichMatrix[insufficientMetaboliteIdx, :] != 0, rxnFluxes > 0))[0]
			rxnNegIdxs = np.where(np.logical_and(self.stoichMatrix[insufficientMetaboliteIdx, :] != 0, rxnFluxes < 0))[0]
			while(np.dot(self.stoichMatrix, rxnFluxes)[insufficientMetaboliteIdx] + moleculeCounts[insufficientMetaboliteIdx] < 0):
				rxnFluxes[rxnPosIdxs] -= 1
				rxnFluxes[rxnNegIdxs] += 1
				rxnFluxes[rxnPosIdxs] = np.fmax(0, rxnFluxes[rxnPosIdxs])
				rxnFluxes[rxnNegIdxs] = np.fmin(0, rxnFluxes[rxnNegIdxs])

		assert(np.all(moleculeCounts + np.dot(self.stoichMatrix, rxnFluxes) >= 0))

		self.molecules.countsInc(
			np.dot(self.stoichMatrix, rxnFluxes)
			)



