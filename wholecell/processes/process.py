"""
Process

Process submodel base class. Defines interface that processes expose to the simulation and to the states.

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/2/2013
"""

from __future__ import absolute_import, division, print_function

import warnings

import wholecell.states.bulk_molecules
import wholecell.states.unique_molecules
import wholecell.states.local_environment
from wholecell.containers.unique_objects_container import Access

import numpy as np

from wholecell.listeners.listener import WriteMethod

class Process(object):
	""" Process """

	_name = 'Process'
	EDIT_ACCESS = (Access.EDIT, )
	EDIT_DELETE_ACCESS = (Access.EDIT, Access.DELETE)

	# Constructor
	def __init__(self):
		# Constants
		self._processIndex = None

		# Simulation reference (used to access time)
		self._sim = None

		# Simulation random stream
		self.randomState = None # set by Simulation
		self.seed = None # set by Simulation

		# References to state
		self._internal_states = None
		self._external_states = None


	# Construct object graph, calculate constants
	def initialize(self, sim, sim_data):
		self._sim = sim

		self._processIndex = sim.processes.keys().index(self._name)

		self._internal_states = sim.internal_states
		self._external_states = sim.external_states


	# Set state partitioning options
	# TODO: make this logic consistent amongst states and allow more options
	def bulkMoleculesRequestPriorityIs(self, priorityLevel):
		self._internal_states["BulkMolecules"].processRequestPriorityIs(
			self._processIndex, priorityLevel)

	def timeStepSec(self):
		return self._sim.timeStepSec()

	def isTimeStepShortEnough(self, *args):
		return True

	def wasTimeStepShortEnough(self, *args):
		return True

	# Construct views
	def environmentView(self, moleculeIDs):
		return wholecell.states.local_environment.EnvironmentView(
			self._external_states['Environment'], self, moleculeIDs)


	def bulkMoleculesView(self, moleculeIDs):
		return wholecell.states.bulk_molecules.BulkMoleculesView(
			self._internal_states['BulkMolecules'], self, moleculeIDs)


	def bulkMoleculeView(self, moleculeIDs):
		return wholecell.states.bulk_molecules.BulkMoleculeView(
			self._internal_states['BulkMolecules'], self, moleculeIDs)


	def uniqueMoleculesView(self, moleculeName):
		return wholecell.states.unique_molecules.UniqueMoleculesView(
			self._internal_states['UniqueMolecules'], self, moleculeName)


	# Communicate with listeners

	# TODO: consider an object-oriented interface to reading/writing to listeners
	# that way, processes would use object handles instead of strings
	def writeToListener(self, listenerName, attributeName, value, writeMethod = WriteMethod.update):
		if listenerName not in self._sim.listeners.viewkeys():
			warnings.warn("The {} process attempted to write {} to the {} listener, but there is no listener with that name.".format(
				self._name,
				attributeName,
				listenerName
				))

		else:
			listener = self._sim.listeners[listenerName]

			if not hasattr(listener, attributeName):
				warnings.warn("The {} process attempted to write {} to the {} listener, but the listener does not have that attribute.".format(
					self._name,
					attributeName,
					listenerName
					))
			else:
				if writeMethod == WriteMethod.update:
					setattr(listener, attributeName, value)
				elif writeMethod == WriteMethod.increment:
					setattr(listener, attributeName, getattr(listener, attributeName) + value)
				elif writeMethod == WriteMethod.append:
					data = getattr(listener, attributeName)
					if isinstance(value, np.ndarray):
						setattr(listener, attributeName, np.append(data, value, axis=0))
					else:
						warnings.warn("The {} process attempted to append to {} on the {} listener, but the given value is not an ndarray".format(
							self._name,
							attributeName,
							listenerName))
				elif writeMethod == WriteMethod.fill:
					data = getattr(listener, attributeName)
					if isinstance(value, np.ndarray) and len(value.shape) == 1:
						n_elements = value.size
						data[:n_elements] = value
						data[n_elements:] = np.nan
					else:
						warnings.warn("The {} process attempted to fill in {} on the {} listener, but the given value is not a 1-dimensional ndarray".format(
							self._name,
							attributeName,
							listenerName))
				else:
					raise warnings.warn("The {} process attempted to write {} to the {} listener, but used an invalid write method: {}".format(
						self._name,
						attributeName,
						listenerName,
						writeMethod))


	def readFromListener(self, listenerName, attributeName):
		if listenerName not in self._sim.listeners.viewkeys():
			raise Exception("The {} process attempted to read {} from the {} listener, but there is no listener with that name.".format(
				self._name,
				attributeName,
				listenerName
				))

		else:
			listener = self._sim.listeners[listenerName]

			if not hasattr(listener, attributeName):
				raise Exception("The {} process attempted to read {} from the {} listener, but the listener does not have that attribute.".format(
					self._name,
					attributeName,
					listenerName
					))

			else:
				return getattr(listener, attributeName)


	# Calculate requests for a single time step
	def calculateRequest(self):
		# Implemented by subclass
		pass


	# Calculate submodel contribution to temporal evolution of cell
	def evolveState(self):
		# Implemented by subclass
		pass


	# Basic accessors

	def time(self):
		return self._sim.time()


	@classmethod
	def name(cls):
		return cls._name
