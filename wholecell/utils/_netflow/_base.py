
from __future__ import division

class NetworkFlowProblemBase(object):
	_maximize = True

	def flowMaterialCoeffIs(self, flow, material, coefficient):
		raise NotImplementedError()

	def flowLowerBoundIs(self, flow, lowerBound):
		raise NotImplementedError()

	def flowUpperBoundIs(self, flow, upperBound):
		raise NotImplementedError()

	def flowObjectiveCoeffIs(self, flow, coefficient):
		raise NotImplementedError()

	def flowRates(self, flows):
		raise NotImplementedError()

	def maximizeObjective(self, doMax):
		self._maximize = doMax
