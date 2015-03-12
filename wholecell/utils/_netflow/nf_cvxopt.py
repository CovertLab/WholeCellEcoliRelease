
from __future__ import division

import numpy as np
import cvxopt
import cvxopt.solvers

from ._base import NetworkFlowProblemBase

class NetworkFlowCvxopt(NetworkFlowProblemBase):
	_upperBoundDefault = np.inf
	_lowerBoundDefault = 0
	_numericalInfinity = 1e6

	def __init__(self):
		self._materialNames = {}
		self._flowNames = {}

		self._materialIndexes = []
		self._flowIndexes = []
		self._coefficients = []

		self._objectiveIndexes = []
		self._objectiveCoefficients = []

		self._lowerBoundIndexes = []
		self._lowerBoundValues = []

		self._upperBoundIndexes = []
		self._upperBoundValues = []

		self._changedEqConst = True
		self._changedIneqConst = True
		self._changedObjective = True
		self._solved = False


	def _getMaterialIndex(self, material):
		try:
			return self._materialNames[material]

		except KeyError:
			self._materialNames[material] = len(self._materialNames)
			return len(self._materialNames) - 1


	def _getFlowIndex(self, flow):
		try:
			return self._flowNames[flow]

		except KeyError:
			self._flowNames[flow] = len(self._flowNames)
			return len(self._flowNames) - 1


	def flowMaterialCoeffIs(self, flow, material, coefficient):
		self._materialIndexes.append(self._getMaterialIndex(material))
		self._flowIndexes.append(self._getFlowIndex(flow))
		self._coefficients.append(coefficient)

		self._changedEqConst = True


	def flowLowerBoundIs(self, flow, lowerBound):
		self._lowerBoundIndexes.append(self._getFlowIndex(flow))
		self._lowerBoundValues.append(lowerBound)

		self._changedIneqConst = True


	def flowUpperBoundIs(self, flow, upperBound):
		self._upperBoundIndexes.append(self._getFlowIndex(flow))
		self._upperBoundValues.append(upperBound)

		self._changedIneqConst = True


	def flowObjectiveCoeffIs(self, flow, coefficient):
		self._objectiveIndexes.append(self._getFlowIndex(flow))
		self._objectiveCoefficients.append(coefficient)

		self._changedObjective = True


	def flowRates(self, flows):
		if isinstance(flows, basestring):
			flows = (flows,)

		self._solve()

		indexes = [self._getFlowIndex(flow) for flow in flows]

		return self._flowRates[indexes]


	def _solve(self):
		self._formEqConst()
		self._formIneqConst()
		self._formObjective()

		if self._solved:
			return

		oldOptions = cvxopt.solvers.options.copy()
		cvxopt.solvers.options["LPX_K_MSGLEV"] = 0

		solution = cvxopt.solvers.lp(self._f, self._G, self._h, self._A, self._b, solver = "glpk")

		cvxopt.solvers.options.update(oldOptions)

		self._rawSolution = solution
		self._flowRates = np.array(self._rawSolution["x"]).flatten()

		self._solved = True


	def _formEqConst(self):
		if not self._changedEqConst:
			return

		self._A = cvxopt.spmatrix(
			self._coefficients,
			self._materialIndexes,
			self._flowIndexes
			)

		self._b = cvxopt.matrix(np.zeros(self._nMaterials(), np.float64))

		self._changedEqConst = False
		self._solved = False


	def _formIneqConst(self):
		if not self._changedIneqConst:
			return

		h = np.empty(2*self._nFlows())

		h[:self._nFlows()] = self._upperBoundDefault
		h[self._nFlows():] = self._lowerBoundDefault

		# TODO: explore whether omitting "infinite" constraints improves
		# solver time/stability

		np.clip(h, -self._numericalInfinity, self._numericalInfinity, h)

		h[self._upperBoundIndexes] = np.array(self._upperBoundValues)
		h[np.array(self._lowerBoundIndexes) + self._nFlows()] = -np.array(self._lowerBoundValues)

		self._h = cvxopt.matrix(h)

		self._G = cvxopt.spmatrix(
			[+1]*self._nFlows() + [-1]*self._nFlows(),
			range(2*self._nFlows()),
			range(self._nFlows()) + range(self._nFlows())
			)

		self._changedIneqConst = False
		self._solved = False


	def _formObjective(self):
		if not self._changedObjective:
			return

		f = np.zeros(self._nFlows())

		f[self._objectiveIndexes] = self._objectiveCoefficients

		self._f = cvxopt.matrix(
			-f if self._maximize else f
			)

		self._changedObjective = False
		self._solved = False


	def _nMaterials(self):
		return len(self._materialNames)


	def _nFlows(self):
		return len(self._flowNames)
