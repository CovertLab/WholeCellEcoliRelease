
from __future__ import division

from collections import defaultdict

from wholecell.utils._netflow import glpk
import numpy as np
from scipy.sparse import coo_matrix

from ._base import NetworkFlowProblemBase

class NetworkFlowGLPK(NetworkFlowProblemBase):
	_lowerBoundDefault = 0
	_upperBoundDefault = np.inf

	def __init__(self):
		self._model = glpk.glpk()
		self._model.set_quiet()

		self._flows = {}
		self._lb = {}
		self._ub = {}
		self._materialCoeffs = defaultdict(list)

		self._eqConstFixed = False
		self._solved = False


	def _getVar(self, flow):
		if flow in self._flows:
			idx = self._flows[flow]
		else:
			self._model.add_cols(1)
			idx = len(self._flows)
			self._lb[flow] = self._lowerBoundDefault
			self._ub[flow] = self._upperBoundDefault
			self._model.set_col_bounds(
				1 + idx,			# GLPK does 1-indexing
				self._lb[flow],
				self._ub[flow],
				)

			self._flows[flow] = idx

		return idx


	def flowMaterialCoeffIs(self, flow, material, coefficient):
		if self._eqConstFixed:
			raise Exception("Cannot add new flow-material coefficients.")

		idx = self._getVar(flow)

		self._materialCoeffs[material].append((coefficient, idx))

		self._solved = False


	def flowLowerBoundIs(self, flow, lowerBound):
		idx = self._getVar(flow)
		self._lb[flow] = lowerBound
		self._model.set_col_bounds(
			1 + idx,				# GLPK does 1 indexing
			self._lb[flow],
			self._ub[flow],
			)

		self._solved = False


	def flowUpperBoundIs(self, flow, upperBound):
		idx = self._getVar(flow)
		self._ub[flow] = upperBound
		self._model.set_col_bounds(
			1 + idx,				# GLPK does 1 indexing
			self._lb[flow],
			self._ub[flow],
			)

		self._solved = False


	def flowObjectiveCoeffIs(self, flow, coefficient):
		idx = self._getVar(flow)
		self._model.set_obj_coef(
			1 + idx,				# GLPK does 1 indexing
			coefficient
			)

		self._solved = False


	def flowRates(self, flows):
		if isinstance(flows, basestring):
			flows = (flows,)

		self._solve()

		return np.array(
			[self._model.get_primal_value(1 + self._getVar(flow)) for flow in flows]
			)


	def _solve(self):
		if self._solved:
			return

		if not self._eqConstFixed:
			self._model.add_rows(len(self._materialCoeffs))
			A = np.zeros((len(self._materialCoeffs), len(self._flows)))
			# avoid creating duplicate constraints
			for materialIdx, (material, pairs) in enumerate(sorted(self._materialCoeffs.viewitems())):
				for pair in pairs:
					A[materialIdx, pair[1]] = pair[0]
			A_coo = coo_matrix(A)
			row = np.hstack((-1, 1 + A_coo.row.astype(np.int32))).astype(np.int32)
			col = np.hstack((-1, 1 + A_coo.col.astype(np.int32))).astype(np.int32)
			data = np.hstack((np.nan, A_coo.data.astype(np.float64))).astype(np.float64)
			self._model.add_eq_constrs(row, col, data)

			self._eqConstFixed = True

		if self._maximize:
			self._model.set_sense_max()
		else:
			self._model.set_sense_min()
		
		self._model.optimize()

		self._solved = True
