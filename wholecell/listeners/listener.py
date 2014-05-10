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
	def update(self):
		return


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


class FakeListener(Listener):
	_name = "FakeListener"

	def __init__(self):
		super(FakeListener, self).__init__()

		print '__init__'


	def initialize(self, sim, kb):
		super(FakeListener, self).initialize(sim, kb)

		print 'initialize'


	def allocate(self):
		super(FakeListener, self).allocate()

		print 'allocate'


	def update(self):
		super(FakeListener, self).update()

		print 'update {}'.format(self.time())
