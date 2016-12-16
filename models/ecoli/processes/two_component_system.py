#!/usr/bin/env python

"""
Two component system

Two component system sub-model

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/3/2016

@author: Heejo Choi
"""

from __future__ import division

import numpy as np
import scipy.integrate

import wholecell.processes.process
from wholecell.utils import units
from wholecell.utils.constants import REQUEST_PRIORITY_TWO_COMPONENT_SYSTEM

import theano.tensor as T
import theano



class TwoComponentSystem(wholecell.processes.process.Process):
	""" Two component system """

	_name = "TwoComponentSystem"

	# Constructor
	def __init__(self):

		super(TwoComponentSystem, self).__init__()


	# Construct object graph
	def initialize(self, sim, sim_data):
		super(TwoComponentSystem, self).initialize(sim, sim_data)

		self.nAvogadro = sim_data.constants.nAvogadro.asNumber(1 / units.mmol)
		self.cellDensity = sim_data.constants.cellDensity.asNumber(units.g / units.L)

		# Create matrices and vectors

		self.stoichMatrix = sim_data.process.two_component_system.stoichMatrix().astype(np.int64)
		self.Rp = sim_data.process.two_component_system.Rp
		self.Pp = sim_data.process.two_component_system.Pp
		self.derivatives = sim_data.process.two_component_system.derivatives
		self.derivativesJacobian = sim_data.process.two_component_system.derivativesJacobian
		self.metsToRxnFluxes = sim_data.process.two_component_system.metsToRxnFluxes
		self.moleculesToNextTimeStep = sim_data.process.two_component_system.moleculesToNextTimeStep

		# Build views

		self.moleculeNames = sim_data.process.two_component_system.moleculeNames
		self.molecules = self.bulkMoleculesView(self.moleculeNames)

		# Set priority to a lower value (but greater priority than metabolism)
		self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_TWO_COMPONENT_SYSTEM)

	def calculateRequest(self):
		moleculeCounts = self.molecules.total()

		cellMass = (self.readFromListener("Mass", "cellMass") * units.fg).asNumber(units.g)
		self.cellVolume = cellMass / self.cellDensity

		self.req, self.allMoleculeChanges = self.moleculesToNextTimeStep(moleculeCounts, self.cellVolume, self.nAvogadro, self.timeStepSec())

		self.molecules.requestIs(self.req)


	def evolveState(self):
		moleculeCounts = self.molecules.counts()

		if (self.req > moleculeCounts).any():
			_, self.allMoleculeChanges = self.moleculesToNextTimeStep(moleculeCounts, self.cellVolume, self.nAvogadro, self.timeStepSec())
			self.molecules.countsInc(self.allMoleculeChanges)
		else:
			self.molecules.countsInc(self.allMoleculeChanges)