"""
Two component system

Two component system sub-model

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/3/2016

"""
from __future__ import absolute_import, division, print_function
import numpy as np

import wholecell.processes.process
from wholecell.utils import units
from wholecell.utils.constants import REQUEST_PRIORITY_TWO_COMPONENT_SYSTEM


class TwoComponentSystem(wholecell.processes.process.Process):
	""" Two component system """

	_name = "TwoComponentSystem"

	# Constructor
	def __init__(self):

		super(TwoComponentSystem, self).__init__()


	# Construct object graph
	def initialize(self, sim, sim_data):
		super(TwoComponentSystem, self).initialize(sim, sim_data)

		# Get constants
		self.nAvogadro = sim_data.constants.nAvogadro.asNumber(1 / units.mmol)
		self.cellDensity = sim_data.constants.cellDensity.asNumber(units.g / units.L)

		# Create method
		self.moleculesToNextTimeStep = sim_data.process.two_component_system.moleculesToNextTimeStep

		# Build views
		self.moleculeNames = sim_data.process.two_component_system.moleculeNames
		self.molecules = self.bulkMoleculesView(self.moleculeNames)

		# Set priority to a lower value (but greater priority than metabolism)
		self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_TWO_COMPONENT_SYSTEM)


	def calculateRequest(self):
		# Get molecule counts
		moleculeCounts = self.molecules.total()

		# Get cell mass and volume
		cellMass = (self.readFromListener("Mass", "cellMass") * units.fg).asNumber(units.g)
		self.cellVolume = cellMass / self.cellDensity

		# Solve ODEs to next time step using the BDF solver through solve_ivp.
		# Note: the BDF solver has been empirically tested to be the fastest
		# solver for this setting among the list of solvers that can be used
		# by the scipy ODE suite.
		self.molecules_required, self.all_molecule_changes = self.moleculesToNextTimeStep(
			moleculeCounts, self.cellVolume, self.nAvogadro,
			self.timeStepSec(), solver="BDF"
			)

		# Request counts of molecules needed
		self.molecules.requestIs(self.molecules_required)


	def evolveState(self):
		# Get counts of molecules allocated to this process
		moleculeCounts = self.molecules.counts()

		# Check if any molecules were allocated fewer counts than requested
		if (self.molecules_required > moleculeCounts).any():

			# Solve ODEs to next time step using the the counts of molecules
			# allocated to this process using the LSODA solver through odeint.
			# Note: for this setting where the counts of molecules are
			# relatively small, the default LSODA solver used by odeint was
			# empirically tested to be the fastest.
			_, self.all_molecule_changes = self.moleculesToNextTimeStep(
				moleculeCounts, self.cellVolume, self.nAvogadro,
				self.timeStepSec()
				)

			# Increment changes in molecule counts
			self.molecules.countsInc(self.all_molecule_changes)
		else:
			# Increment changes in molecule counts
			self.molecules.countsInc(self.all_molecule_changes)
