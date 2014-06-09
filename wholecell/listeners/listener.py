"""
Listener

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/10/2014
"""

class Listener(object):
	_name = None

	def __init__(self):
		pass

	# Construct state-process graph, calculate constants
	def initialize(self, sim, kb):
		self._sim = sim

		self._nProcesses = len(sim.processes)


	# Allocate memory
	def allocate(self):
		pass


	# Calculate (and cache) any dependent properties
	def initialUpdate(self):
		# Default behavior is to call the normal update method
		
		self.update()


	def updatePostRequest(self):
		pass


	def update(self):
		pass


	# Saving and loading

	def pytablesCreate(self, h5file, expectedRows):
		pass

	def pytablesAppend(self, h5file):
		pass

	def pytablesLoad(self, h5file, timePoint):
		pass


	# Basic accessors

	def time(self):
		return self._sim.time()


	def timeStep(self):
		return self._sim.timeStep()


	@classmethod
	def name(cls):
		return cls._name
