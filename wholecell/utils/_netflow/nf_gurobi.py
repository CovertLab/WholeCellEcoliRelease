
from __future__ import division

from collections import defaultdict

import gurobipy as grb
import numpy as np

from ._base import NetworkFlowProblemBase

class NetworkFlowGurobi(NetworkFlowProblemBase):
	_lowerBoundDefault = 0
	_upperBoundDefault = grb.GRB.INFINITY

	def __init__(self):
		self._model = grb.Model()
		self._model.setParam("OutputFlag", False)

		self._flows = {}
		self._materialCoeffs = defaultdict(list)
		self._objectiveCoefficients = []

		self._eqConstFixed = False
		self._solved = False


	def _getVar(self, flow):
		try:
			var = self._flows[flow]

		except KeyError:
			var = self._model.addVar(
				name = flow,
				lb = self._lowerBoundDefault,
				ub = self._upperBoundDefault
				)

			self._model.update()

			self._flows[flow] = var

		return var


	def flowMaterialCoeffIs(self, flow, material, coefficient):
		if self._eqConstFixed:
			raise Exception("Cannot add new flow-material coefficients.")

		var = self._getVar(flow)

		self._materialCoeffs[material].append((coefficient, var))

		self._solved = False


	def flowLowerBoundIs(self, flow, lowerBound):
		var = self._getVar(flow)
		var.setAttr("lb", lowerBound)

		self._solved = False


	def flowUpperBoundIs(self, flow, upperBound):
		var = self._getVar(flow)
		var.setAttr("ub", upperBound)

		self._solved = False


	def flowObjectiveCoeffIs(self, flow, coefficient):
		var = self._getVar(flow)
		self._objectiveCoefficients.append((coefficient, var))

		self._solved = False


	def flowRates(self, flows):
		if isinstance(flows, basestring):
			flows = (flows,)

		self._solve()

		return np.array(
			[self._getVar(flow).X for flow in flows]
			)


	def _solve(self):
		if self._solved:
			return

		if not self._eqConstFixed:
			# avoid creating duplicate constraints
			for material, pairs in self._materialCoeffs.viewitems():
				self._model.addConstr(
					grb.LinExpr(pairs) == 0,
					material
					)

			self._eqConstFixed = True

		self._model.setObjective(
			grb.LinExpr(self._objectiveCoefficients),
			grb.GRB.MAXIMIZE if self._maximize else grb.GRB.MINIMIZE
			)

		self._model.update()

		self._model.optimize()

		self._solved = True
