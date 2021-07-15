# noinspection HttpUrlsUsage
"""
NetworkFlow's interface to GNU Linear Programming Kit.
This uses a fraction of the API that's available via the swiglpk pip.

The GLPK manual is at http://kam.mff.cuni.cz/~elias/glpk.pdf

swiglpk is a low-level, SWIG-generated interface that mechanically
bridges the gap between Python and GLPK's C API.
See https://pypi.org/project/swiglpk/ and
https://github.com/biosustain/swiglpk

See http://www.swig.org/Doc1.3/Python.html for documentation on SWIG use
from Python. E.g. an intArray is just a memory pointer that can get and
set native int values. It doesn't even track its length, do bounds
checks, or support iteration, so be careful!

CAUTION: GLPK will exit the process on error, so do pre-checks. It's
possible to set its error_hook -- if SWIG can pass in a native function
pointer and use setjmp/longjmp, and the hook has to reallocate the GLPK
environment, losing state info.
"""

# TODO(Jerry): Check that integer arguments are in range so GLPK won't exit?

from __future__ import absolute_import, division, print_function

from collections import defaultdict
from enum import Enum

import numpy as np
from scipy.sparse import coo_matrix
import six
import swiglpk as glp

from ._base import NetworkFlowProblemBase
from six.moves import range
from six.moves import zip

class MessageLevel(Enum):
	OFF = glp.GLP_MSG_OFF  # no output
	ERR = glp.GLP_MSG_ERR  # error and warning messages only (default)
	ON  = glp.GLP_MSG_ON   # normal output
	ALL = glp.GLP_MSG_ALL  # full output (including informational messages)
	DBG = glp.GLP_MSG_DBG

class SimplexMethod(Enum):
	PRIMAL = glp.GLP_PRIMAL  # two-phase primal simplex (the default)
	DUALP = glp.GLP_DUALP    # two-phase dual simplex
	DUAL = glp.GLP_DUAL      # two-phase dual simplex; fallback to the primal simplex

SOLUTION_STATUS_TO_STRING = {
    glp.GLP_OPT:	'GLP_OPT: optimal',  # solution is optimal
    glp.GLP_FEAS:	'GLP_FEAS: feasible',  # solution is feasible
    glp.GLP_INFEAS:	'GLP_INFEAS: infeasible',  # solution is infeasible
    glp.GLP_NOFEAS:	'GLP_NOFEAS: no feasible',  # problem has no feasible solution
    glp.GLP_UNBND:	'GLP_UNBND: unbounded',  # problem has no unbounded solution
    glp.GLP_UNDEF:	'GLP_UNDEF: undefined',  # solution is undefined
}

_MAXED_OUT = ("GLP_EOBJ{}L: Dual simplex: The objective function being"
			  + " maximized reached its {} limit and continues {}")

SIMPLEX_RETURN_CODE_TO_STRING = {
	0: '',  # successfully solved; not necessarily an optimal solution
	glp.GLP_EBADB:  "GLP_EBADB: Basis is invalid, number of basic variables != number of rows",
	glp.GLP_ESING:  "GLP_ESING: Basis matrix is singular",
	glp.GLP_ECOND:  "GLP_ECOND: Basis matrix is ill-conditioned, it's condition number is too large",
	glp.GLP_EBOUND: "GLP_EBOUND: Some double-bounded variables have incorrect bounds",
	glp.GLP_EFAIL:  "GLP_EFAIL: Solver failure",
	glp.GLP_EOBJLL: _MAXED_OUT.format('L', 'lower', 'decreasing'),
	glp.GLP_EOBJUL: _MAXED_OUT.format('U', 'upper', 'increasing'),
	glp.GLP_EITLIM: "GLP_EITLIM: Iteration limit exceeded",
	glp.GLP_ETMLIM: "GLP_ETMLIM: Time limit exceeded",
	glp.GLP_ENOPFS: "GLP_ENOPFS: Presolver: Problem has no primal feasible solution",
	glp.GLP_ENODFS: "GLP_ENODFS: Presolver: Problem has no dual feasible solution",
}

def _toIndexArray(array):
	"""Convert an array to a GLPK IntArray of indexes: Convert the indexes to
	int, add 1 to each, and prepend a dummy value.
	"""

	ia = glp.intArray(len(array) + 1)
	ia[0] = -1
	for (i, value) in enumerate(array):
		ia[i + 1] = int(value) + 1
	return ia

def _toDoubleArray(array):
	"""Convert an array to a GLPK DoubleArray of indexes: Convert the values to
	double and prepend a dummy value.
	"""

	da = glp.doubleArray(len(array) + 1)
	da[0] = np.nan
	for (i, value) in enumerate(array):
		da[i + 1] = float(value)
	return da


class NetworkFlowGLPK(NetworkFlowProblemBase):
	def __init__(self, quadratic_objective=False):
		if quadratic_objective:
			raise ValueError('Quadratic objective not supported for GLPK')

		self._lp = glp.glp_create_prob()
		self._smcp = glp.glp_smcp()  # simplex solver control parameters
		glp.glp_init_smcp(self._smcp)
		self._smcp.msg_lev = glp.GLP_MSG_ERR
		self.simplex_iteration_limit = 10000
		self._n_vars = 0
		self._n_eq_constraints = 0

		self._flows = {}
		self._lb = {}
		self._ub = {}
		self._objective = {}
		self._materialCoeffs = defaultdict(list)
		self._materialIdxLookup = {}
		self._flow_index_arrays = {}
		self._coeff_arrays = {}
		self._flow_locations = {}

		self._eqConstBuilt = False
		self._solved = False

		self.inf = np.inf

		self._lowerBoundDefault = 0
		self._upperBoundDefault = self.inf

	def __del__(self):
		glp.glp_delete_prob(self._lp)


	@property
	def message_level(self):
		"""The message level for terminal output, as an enum value."""
		return MessageLevel(self._smcp.msg_lev)

	@message_level.setter
	def message_level(self, message_level):
		"""Set the message level from a MessageLevel enum value."""
		self._smcp.msg_lev = message_level.value

	@property
	def simplex_method(self):
		"""The Simplex method option, as an enum value."""
		return SimplexMethod(self._smcp.meth)

	@simplex_method.setter
	def simplex_method(self, simplex_method):
		"""Set the Simplex method option from a SimplexMethod enum value."""
		self._smcp.meth = simplex_method.value

	@property
	def simplex_iteration_limit(self):
		"""The Simplex iteration limit, an integer."""
		return self._smcp.it_lim

	@simplex_iteration_limit.setter
	def simplex_iteration_limit(self, limit):
		"""Set the Simplex iteration limit."""
		self._smcp.it_lim = int(limit)

	@property
	def primal_feasible_tolerance(self):
		"""Tolerance used to check if the basic solution is primal feasible."""
		return self._smcp.tol_bnd

	@primal_feasible_tolerance.setter
	def primal_feasible_tolerance(self, tolerance):
		"""Tolerance used to check if the basic solution is primal feasible.
		(Do not change this parameter without detailed understanding its purpose.)
		Default: 1e-7.
		"""
		self._smcp.tol_bnd = float(tolerance)

	@property
	def status_code(self):
		"""The generic status code for the current basic solution."""
		return glp.glp_get_status(self._lp)

	@property
	def status_string(self):
		"""Return the generic status message for the current basic solution."""
		return SOLUTION_STATUS_TO_STRING.get(
			self.status_code, "GLP_?: UNKNOWN SOLUTION STATUS CODE")

	def columnPrimalValue(self, j):
		"""Return the primal value of the structural variable for j-th column."""
		return glp.glp_get_col_prim(self._lp, j)

	def columnDualValue(self, j):
		"""Return the dual value (i.e. reduced cost) of the structural variable
		for the j-th column.
		"""
		return glp.glp_get_col_dual(self._lp, j)

	def rowPrimalValue(self, i):
		"""Return the primal value of the structural variable for i-th row."""
		return glp.glp_get_row_prim(self._lp, i)

	def rowDualValue(self, i):
		"""Return the dual value (i.e. reduced cost) of the structural variable
		for the i-th row.
		"""
		return glp.glp_get_row_dual(self._lp, i)


	def _add_rows(self, n_rows):
		"""Add n_rows rows (constraints) to the end of the row list. Each new
		row is initially free (unbounded) and has an empty list of constraint
		coefficients.
		"""
		glp.glp_add_rows(self._lp, n_rows)
		self._n_eq_constraints += n_rows

	def _add_cols(self, n_cols):
		"""Add n_cols columns (structural variables) to the end of the column
		list. Each new column is initially fixed at zero and has an empty list
		of constraint coefficients.
		"""
		glp.glp_add_cols(self._lp, n_cols)
		self._n_vars += n_cols

	def _set_col_bounds(self, index, lower, upper):
		"""Set the type and bounds of index-th (j-th) column (structural
		variable). Use -np.inf or np.inf to specify no lower or upper bound.
		"""
		if np.isinf(lower) and np.isinf(upper):
			variable_type = glp.GLP_FR  # free (unbounded) variable
		elif lower == upper:
			variable_type = glp.GLP_FX  # fixed variable
		elif lower > upper:
			raise ValueError("The lower bound must be <= upper bound")
		elif not np.isinf(lower) and not np.isinf(upper):
			variable_type = glp.GLP_DB  # double-bounded variable
		elif np.isinf(upper):
			variable_type = glp.GLP_LO  # variable with lower bound
		else:
			variable_type = glp.GLP_UP  # variable with upper bound

		glp.glp_set_col_bnds(self._lp, index, variable_type, lower, upper)

	def _getVar(self, flow):
		if flow in self._flows:
			idx = self._flows[flow]
		elif self._eqConstBuilt:
			raise ValueError('Equality constraints already built. Unable to add new flow: "{}".'.format(flow))
		else:
			self._add_cols(1)
			idx = len(self._flows)
			self._lb[flow] = self._lowerBoundDefault
			self._ub[flow] = self._upperBoundDefault
			self._set_col_bounds(
				1 + idx,			# GLPK does 1-indexing
				self._lb[flow],
				self._ub[flow],
				)

			self._flows[flow] = idx

		return idx

	def setFlowMaterialCoeff(self, flow, material, coefficient):
		if self._eqConstBuilt:
			if material not in self._materialIdxLookup:
				raise ValueError("Invalid material: {}".format(material))
			if flow not in self._flows:
				raise ValueError("Invalid flow: {}".format(flow))

			length = len(self._flow_locations[material])
			flow_loc = self._flow_locations[material][self._flows[flow]]
			rowIdx = int(self._materialIdxLookup[material] + 1)
			colIdxs = self._flow_index_arrays[material]
			data = self._coeff_arrays[material]
			data[flow_loc] = float(coefficient)  # swiglpk offsets index by 1
			glp.glp_set_mat_row(self._lp, rowIdx, length, colIdxs, data)
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

		reset_bounds = False
		idx = self._getVar(flow)
		if lowerBound is not None and self._lb[flow] != lowerBound:
			self._lb[flow] = lowerBound
			# Need to prevent cases where the new lower bound is greater than
			# current upper bound and a new upper bound is not provided here
			if upperBound is None and self._lb[flow] > self._ub[flow]:
				self._ub[flow] = self._lb[flow]
			reset_bounds = True
		if upperBound is not None and self._ub[flow] != upperBound:
			self._ub[flow] = upperBound
			# Need to prevent cases where the new upper bound is less than
			# current lower bound and a new loewr bound is not provided here
			if lowerBound is None and self._lb[flow] > self._ub[flow]:
				self._lb[flow] = self._ub[flow]
			reset_bounds = True

		if reset_bounds:
			self._set_col_bounds(
				1 + idx,  # GLPK does 1 indexing
				self._lb[flow],
				self._ub[flow],
				)
			self._solved = False

	def setFlowObjectiveCoeff(self, flow, coefficient):
		idx = self._getVar(flow)
		self._objective[flow] = coefficient
		glp.glp_set_obj_coef(
			self._lp,
			1 + idx,  # GLPK does 1 indexing
			coefficient
			)

		self._solved = False

	def getFlowObjectiveCoeff(self, flow):
		return self._objective[flow]

	def getFlowRates(self, flows):
		if isinstance(flows, six.string_types):
			flows = (flows,)

		self._solve()

		return np.array(
			[self._col_primals[self._flows[flow]]
			if flow in self._flows else None
			for flow in flows]
			)

	def getShadowPrices(self, materials):
		if not self._eqConstBuilt:
			raise RuntimeError("Equality constraints not yet built. Finish construction of the problem before accessing dual values.")

		self._solve()

		return np.array(
			[self._row_duals[self._materialIdxLookup[material]]
			if material in self._materialIdxLookup else None
			for material in materials]
			)

	def getReducedCosts(self, fluxNames):
		if not self._eqConstBuilt:
			raise RuntimeError("Equality constraints not yet built. Finish construction of the problem before accessing dual values.")

		self._solve()

		return np.array(
			[self._col_duals[self._flows[fluxName]]
			if fluxName in self._flows else None
			for fluxName in fluxNames]
			)

	def getObjectiveValue(self):
		"""The current value of the objective function."""
		self._solve()
		return glp.glp_get_obj_val(self._lp)

	def getSMatrix(self):
		if not self._eqConstBuilt:
			raise RuntimeError("Equality constraints not yet built. Finish construction of the problem before accessing S matrix.")
		A = np.zeros((len(self._materialCoeffs), len(self._flows)))
		for materialIdx, material in enumerate(sorted(self._materialCoeffs)):
			coeffs = self._coeff_arrays[material]
			for flow_idx, coeff_idx in self._flow_locations[material].items():
				A[materialIdx, flow_idx] = coeffs[coeff_idx]
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
		n_coeffs = len(self._materialCoeffs)
		n_flows = len(self._flows)

		self._add_rows(n_coeffs)
		A = np.zeros((n_coeffs, n_flows))
		# avoid creating duplicate constraints
		self._materialIdxLookup = {}
		for materialIdx, (material, pairs) in enumerate(sorted(six.viewitems(self._materialCoeffs))):
			self._materialIdxLookup[material] = materialIdx
			for pair in pairs:
				A[materialIdx, pair[1]] = pair[0]

		A_coo = coo_matrix(A)
		rowIdxs = _toIndexArray(A_coo.row)
		colIdxs = _toIndexArray(A_coo.col)
		data = _toDoubleArray(A_coo.data)
		n_elems = len(A_coo.row)

		for row in range(1, self._n_eq_constraints + 1):
			glp.glp_set_row_bnds(self._lp, row, glp.GLP_FX, 0.0, 0.0)
		glp.glp_load_matrix(self._lp, n_elems, rowIdxs, colIdxs, data)

		self._cache_glp_arrays()

		self._eqConstBuilt = True

	def _cache_glp_arrays(self):
		'''
		Caches the arrays needed for the swiglpk interface so they don't need
		to be recomputed at each time step.  Called after the constraints are
		built so the indices won't change.  Coefficients will still need to be
		updated.

		Note: self._flow_index_arrays and self._coeff_arrays mirror the
		information in self._materialCoeffs but also contain padded values for
		the solver. self._materialCoeffs is used in	getSMatrix() and is useful
		for efficient indexing when updating coefficients so it should be kept.
		'''

		for material in self._materialCoeffs:
			coeff, flowIdxs = list(zip(*self._materialCoeffs[material]))
			self._flow_locations[material] = {idx: i + 1  # +1 for swiglpk indexing
				for i, idx in enumerate(flowIdxs)}
			flowIdxs = _toIndexArray(flowIdxs)
			coeff = _toDoubleArray(coeff)
			self._flow_index_arrays[material] = flowIdxs
			self._coeff_arrays[material] = coeff

	def _solve(self):
		if self._solved:
			return

		if self._maximize:
			glp.glp_set_obj_dir(self._lp, glp.GLP_MAX)
		else:
			glp.glp_set_obj_dir(self._lp, glp.GLP_MIN)

		result = glp.glp_simplex(self._lp, self._smcp)

		# Adjust solver options for robustness
		## If no solution within iteration limit, switch to dual method to find solution
		if result == glp.GLP_EITLIM:
			print('Warning: could not find solution with primal method, switching to dual')
			self.simplex_method = SimplexMethod.DUALP
			result = glp.glp_simplex(self._lp, self._smcp)
			self.simplex_method = SimplexMethod.PRIMAL

			## If still no solution, presolve the problem to reduce complexity
			## (also prevents warm start)
			if self.status_code != glp.GLP_OPT:
				print('Warning: could not find solution with dual method, using primal with presolve')
				self._smcp.presolve = glp.GLP_ON
				result = glp.glp_simplex(self._lp, self._smcp)
				self._smcp.presolve = glp.GLP_OFF

		if result != 0:
			raise RuntimeError(SIMPLEX_RETURN_CODE_TO_STRING.get(
				result, "GLP_?: UNKNOWN SOLVER RETURN VALUE"))
		if self.status_code != glp.GLP_OPT:
			raise RuntimeError(self.status_string)

		self._solved = True

		# Read results for better performance when accessing individual values
		self._col_primals = glp.get_col_primals(self._lp)
		self._col_duals = glp.get_col_duals(self._lp)
		self._row_duals = glp.get_row_duals(self._lp)
