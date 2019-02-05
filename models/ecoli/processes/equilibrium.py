#!/usr/bin/env python

"""
Equilibrium

Equilibrium binding sub-model

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 8/21/2015
"""

import numpy as np

from wholecell.utils import units
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

		# Get constants
		self.nAvogadro = sim_data.constants.nAvogadro.asNumber(1 / units.mol)
		self.cellDensity = sim_data.constants.cellDensity.asNumber(units.g / units.L)

		# Create matrix and method
		self.stoichMatrix = sim_data.process.equilibrium.stoichMatrix().astype(np.int64)
		self.fluxesAndMoleculesToSS = sim_data.process.equilibrium.fluxesAndMoleculesToSS
		self.product_indices = [idx for idx in np.where(np.any(self.stoichMatrix > 0, axis=1))[0]]

		# Build views
		moleculeNames = sim_data.process.equilibrium.moleculeNames
		self.molecules = self.bulkMoleculesView(moleculeNames)


	def calculateRequest(self):
		# Get molecule counts
		moleculeCounts = self.molecules.total_counts()

		# Get cell mass and volume
		cellMass = (self.readFromListener("Mass", "cellMass") * units.fg).asNumber(units.g)
		cellVolume = cellMass / self.cellDensity

		# Solve ODEs to steady state
		self.rxnFluxes, self.req = self.fluxesAndMoleculesToSS(moleculeCounts, cellVolume, self.nAvogadro)

		# Request counts of molecules needed
		self.molecules.requestIs(self.req)


	def evolveState(self):
		# Get counts of molecules allocated to this process
		moleculeCounts = self.molecules.counts()

		# If we didn't get allocated all the molecules we need, make do with what we have
		# (decrease reaction fluxes so that they make use of what we have, but not more)
		rxnFluxes = self.rxnFluxes.copy()
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

		# Increment changes in molecule counts
		deltaMolecules = np.dot(self.stoichMatrix, rxnFluxes)
		self.molecules.countsInc(deltaMolecules)

		# Write outputs to listeners
		self.writeToListener("EquilibriumListener", "reactionRates", (
			deltaMolecules[self.product_indices] / self.timeStepSec()
			))
