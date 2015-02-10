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

		if self._sim.loggers.has_key("Shell"):
			self._shellLogger = self._sim.loggers["Shell"]

		else:
			self._shellLogger = None


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

	def tableCreate(self, tableWriter):
		pass

	def tableAppend(self, tableWriter):
		pass

	def tableLoad(self, tableReader, tableIndex):
		pass


	# Basic accessors

	def time(self):
		return self._sim.time()


	def timeStep(self):
		return self._sim.timeStep()


	# Features for shell logging

	def registerLoggedQuantity(self, header, attribute, format_spec = "f", cell_size = 0):
		if self._shellLogger is None:
			return

		self._shellLogger.columnSpecs.append(dict(
			header = header,
			target = "Listener:{}".format(self._name),
			property = attribute,
			format = format_spec,
			length = cell_size,
			sum = False # TODO: get rid of 'sum'
			))


	@classmethod
	def name(cls):
		return cls._name
