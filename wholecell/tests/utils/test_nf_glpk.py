# TODO(Jerry): Fill out test_nf_glpk() and add a performance test.

from __future__ import absolute_import
from __future__ import division

import unittest

import swiglpk as glp

from wholecell.utils._netflow.nf_glpk import (MessageLevel, NetworkFlowGLPK,
	SimplexMethod)


class Test_NetworkFlowGLPK(unittest.TestCase):

	def test_swiglpk(self):
		"""Test the underlying GLPK lib and its SWIG interface based on
		the example from https://github.com/biosustain/swiglpk
		"""
		ia = glp.intArray(1 + 1000)
		ja = glp.intArray(1 + 1000)
		ar = glp.doubleArray(1 + 1000)

		lp = glp.glp_create_prob()
		smcp = glp.glp_smcp()
		glp.glp_init_smcp(smcp)
		smcp.msg_lev = glp.GLP_MSG_ALL  # use GLP_MSG_ERR?

		glp.glp_set_prob_name(lp, "sample")
		glp.glp_set_obj_dir(lp, glp.GLP_MAX)

		glp.glp_add_rows(lp, 3)
		glp.glp_set_row_name(lp, 1, "p")
		glp.glp_set_row_bnds(lp, 1, glp.GLP_UP, 0.0, 100.0)
		glp.glp_set_row_name(lp, 2, "q")
		glp.glp_set_row_bnds(lp, 2, glp.GLP_UP, 0.0, 600.0)
		glp.glp_set_row_name(lp, 3, "r")
		glp.glp_set_row_bnds(lp, 3, glp.GLP_UP, 0.0, 300.0)
		glp.glp_add_cols(lp, 3)
		glp.glp_set_col_name(lp, 1, "x1")
		glp.glp_set_col_bnds(lp, 1, glp.GLP_LO, 0.0, 0.0)
		glp.glp_set_obj_coef(lp, 1, 10.0)
		glp.glp_set_col_name(lp, 2, "x2")
		glp.glp_set_col_bnds(lp, 2, glp.GLP_LO, 0.0, 0.0)
		glp.glp_set_obj_coef(lp, 2, 6.0)
		glp.glp_set_col_name(lp, 3, "x3")
		glp.glp_set_col_bnds(lp, 3, glp.GLP_LO, 0.0, 0.0)
		glp.glp_set_obj_coef(lp, 3, 4.0)

		ia[1] = 1; ja[1] = 1; ar[1] = 1.0  # a[1,1] = 1
		ia[2] = 1; ja[2] = 2; ar[2] = 1.0  # a[1,2] = 1
		ia[3] = 1; ja[3] = 3; ar[3] = 1.0  # a[1,3] = 1
		ia[4] = 2; ja[4] = 1; ar[4] = 10.0  # a[2,1] = 10
		ia[5] = 3; ja[5] = 1; ar[5] = 2.0  # a[3,1] = 2
		ia[6] = 2; ja[6] = 2; ar[6] = 4.0  # a[2,2] = 4
		ia[7] = 3; ja[7] = 2; ar[7] = 2.0  # a[3,2] = 2
		ia[8] = 2; ja[8] = 3; ar[8] = 5.0  # a[2,3] = 5
		ia[9] = 3; ja[9] = 3; ar[9] = 6.0  # a[3,3] = 6

		glp.glp_load_matrix(lp, 9, ia, ja, ar)
		glp.glp_simplex(lp, smcp)

		Z = glp.glp_get_obj_val(lp)
		x1 = glp.glp_get_col_prim(lp, 1)
		x2 = glp.glp_get_col_prim(lp, 2)
		x3 = glp.glp_get_col_prim(lp, 3)

		self.assertAlmostEqual(Z, 733.3333, 4)
		self.assertAlmostEqual(x1, 33.3333, 4)
		self.assertAlmostEqual(x2, 66.6667, 4)
		self.assertAlmostEqual(x3, 0)

		glp.glp_delete_prob(lp)

	def test_nf_glpk(self):
		"""Test the NetworkFlowGLPK interface to GLPK."""
		nf = NetworkFlowGLPK()

		nf.message_level = MessageLevel.ALL
		self.assertEqual(nf.message_level, MessageLevel.ALL)

		self.assertEqual(nf.simplex_method, SimplexMethod.PRIMAL)
		nf.simplex_method = SimplexMethod.DUAL
		self.assertEqual(nf.simplex_method, SimplexMethod.DUAL)

		limit = 1000
		nf.simplex_iteration_limit = limit
		self.assertEqual(nf.simplex_iteration_limit, limit)

		self.assertEqual(nf.status_code, glp.GLP_UNDEF)
		self.assertIn('UNDEF', nf.status_string)

		# TODO: How to use the NetworkFlowGLPK interface?
