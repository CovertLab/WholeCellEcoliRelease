
from __future__ import division

class NetworkFlowProblemBase(object):
	_maximize = True

	def flowMaterialCoeffIs(self, flow, material, coefficient):
		raise NotImplementedError()

	def setFlowBounds(selfs, flow, ub=None, lb=None):
		raise NotImplementedError()

	def flowObjectiveCoeffIs(self, flow, coefficient):
		raise NotImplementedError()

	def flowRates(self, flows):
		raise NotImplementedError()

	def maximizeObjective(self, doMax):
		self._maximize = doMax
