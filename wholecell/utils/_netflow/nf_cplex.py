'''
NetworkFlow interface to the CPLEX solver package. This package supports both
linear and quadratic objectives.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/14/2018
'''

from __future__ import absolute_import, division, print_function

from collections import defaultdict

# NOTE: This file assumes callers catch the ImportError if IBM CPLEX is not
# installed. To use it, install the CPLEX binary library from IBM (it's
# free for students) and do `pip install cplex>=12.8.0.0`.
# noinspection PyPackageRequirements
import cplex
import numpy as np
from scipy.sparse import coo_matrix
import six

from ._base import NetworkFlowProblemBase
from six.moves import zip


class NetworkFlowCPLEX(NetworkFlowProblemBase):
	def __init__(self, quadratic_objective=False):
		self._model = cplex.Cplex()

		self._model.set_log_stream(None)
		self._model.set_error_stream(None)
		self._model.set_warning_stream(None)
		self._model.set_results_stream(None)

		self._flows = {}
		self._lb = {}
		self._ub = {}
		self._objective = {}
		self._materialCoeffs = defaultdict(list)

		self.quadratic_objective = quadratic_objective

		self._eqConstBuilt = False
		self._solved = False

		self.inf = cplex.infinity

		self._lowerBoundDefault = 0
		self._upperBoundDefault = self.inf

	def _getVar(self, flow):
		if flow in self._flows:
			idx = self._flows[flow]
		elif self._eqConstBuilt:
			raise ValueError('Equality constraints already built. Unable to add new flow: "{}".'.format(flow))
		else:
			self._model.variables.add(obj=[0])
			idx = len(self._flows)
			self._lb[flow] = self._lowerBoundDefault
			self._ub[flow] = self._upperBoundDefault
			self._model.variables.set_lower_bounds(idx, self._lb[flow])
			self._model.variables.set_upper_bounds(idx, self._ub[flow])
			self._flows[flow] = idx

		return idx

	def setFlowMaterialCoeff(self, flow, material, coefficient):
		if self._eqConstBuilt:
			materialIdx = self._materialIdxLookup.get(material, None)
			if materialIdx is None:
				raise ValueError("Invalid material: {}".format(material))

			flowIdx = self._flows.get(flow, None)
			if flowIdx is None:
				raise ValueError("Invalid flow: {}".format(flow))

			coeffs, flowIdxs = list(zip(*self._materialCoeffs[material]))
			coeffs = list(coeffs)
			flowLoc = flowIdxs.index(flowIdx)
			coeffs[flowLoc] = coefficient
			self._materialCoeffs[material] = list(zip(coeffs, flowIdxs))

			self._model.linear_constraints.set_coefficients(list(zip([materialIdx], [flowIdx], [coefficient])))
		else:
			idx = self._getVar(flow)
			self._materialCoeffs[material].append((coefficient, idx))

		self._solved = False

	def setFlowBounds(self, flow, lowerBound=None, upperBound=None):
		"""
		Set the lower and upper bounds for a given flow
		inputs:
			flow (str) - name of flow to set bounds for
			lowerBound (float) - lower bound for flow (None if unchanged)
			upperBound (float) - upper bound for flow (None if unchanged)
		"""

		idx = self._getVar(flow)
		if lowerBound is not None:
			self._lb[flow] = lowerBound
		if upperBound is not None:
			self._ub[flow] = upperBound

		self._model.variables.set_lower_bounds(idx, self._lb[flow])
		self._model.variables.set_upper_bounds(idx, self._ub[flow])

		self._solved = False

	def setFlowObjectiveCoeff(self, flow, coefficient):
		idx = self._getVar(flow)
		self._objective[flow] = coefficient

		if self.quadratic_objective:
			self._model.objective.set_quadratic_coefficients(
				idx,
				idx,
				coefficient
				)
		else:
			self._model.objective.set_linear(
				idx,
				coefficient
				)

		self._solved = False

	def getFlowObjectiveCoeff(self, flow):
		return self._objective[flow]

	def getFlowRates(self, flows):
		if isinstance(flows, six.string_types):
			flows = (flows,)

		self._solve()

		flows = [self._flows.get(flow, None) for flow in flows]
		return np.array(self._model.solution.get_values(flows))

	def getShadowPrices(self, materials):
		# reduced cost of row
		if not self._eqConstBuilt:
			raise RuntimeError("Equality constraints not yet built. Finish construction of the problem before accessing dual values.")

		self._solve()

		materials = [self._materialIdxLookup.get(material, None) for material in materials]
		return np.array(self._model.solution.get_dual_values(materials))

	def getReducedCosts(self, fluxNames):
		# reduced cost of col
		if not self._eqConstBuilt:
			raise RuntimeError("Equality constraints not yet built. Finish construction of the problem before accessing dual values.")

		self._solve()

		flows = [self._flows.get(flow, None) for flow in fluxNames]
		return np.array(self._model.solution.get_reduced_costs(flows))

	def getObjectiveValue(self):
		self._solve()
		return self._model.solution.get_objective_value()

	def getSMatrix(self):
		if not self._eqConstBuilt:
			raise RuntimeError("Equality constraints not yet built. Finish construction of the problem before accessing S matrix.")
		A = np.zeros((len(self._materialCoeffs), len(self._flows)))
		self._materialIdxLookup = {}
		for materialIdx, (material, pairs) in enumerate(sorted(six.viewitems(self._materialCoeffs))):
			self._materialIdxLookup[material] = materialIdx
			for pair in pairs:
				A[materialIdx, pair[1]] = pair[0]
		return A

	def getFlowNames(self):
		if not self._eqConstBuilt:
			raise RuntimeError("Equality constraints not yet built. Finish construction of the problem before accessing flow names.")
		return sorted(self._flows, key=self._flows.__getitem__)

	def getMaterialNames(self):
		if not self._eqConstBuilt:
			raise RuntimeError("Equality constraints not yet built. Finish construction of the problem before accessing material names.")
		return sorted(self._materialIdxLookup, key=self._materialIdxLookup.__getitem__)

	def getUpperBounds(self):
		return self._ub.copy()

	def getLowerBounds(self):
		return self._lb.copy()

	def getObjective(self):
		return self._objective.copy()

	def buildEqConst(self):
		if self._eqConstBuilt:
			raise RuntimeError("Equality constraints already built.")

		nMaterials = len(self._materialCoeffs)
		nFlows = len(self._flows)
		A = np.zeros((nMaterials, nFlows))

		# avoid creating duplicate constraints
		self._materialIdxLookup = {}
		for materialIdx, (material, pairs) in enumerate(sorted(six.viewitems(self._materialCoeffs))):
			self._materialIdxLookup[material] = materialIdx
			for pair in pairs:
				A[materialIdx, pair[1]] = pair[0]
		A_coo = coo_matrix(A)
		row = A_coo.row.astype(np.int32).astype(np.int32)
		col = A_coo.col.astype(np.int32).astype(np.int32)
		data = A_coo.data.astype(np.float64).astype(np.float64)

		# set solver constraints
		self._model.linear_constraints.add(rhs=np.zeros(nMaterials), senses='E'*nMaterials)
		self._model.linear_constraints.set_coefficients(list(zip(row.tolist(), col.tolist(), data.tolist())))

		self._eqConstBuilt = True

	def _solve(self):
		if self._solved:
			return

		if self._maximize:
			self._model.objective.set_sense(self._model.objective.sense.maximize)
		else:
			self._model.objective.set_sense(self._model.objective.sense.minimize)

		self._model.solve()

		self._solved = True
