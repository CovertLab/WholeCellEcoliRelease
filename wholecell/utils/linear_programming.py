#!/usr/bin/env python

"""
linearProgramming
Wrapper over cvxopt to call glpk (can be extended in the future to use other solvers)
Solves the problem:

maximize f' x
subject to
	Ax == b
	lb <= x <= ub

Example call (trivial solution):

>>> wholecell.util.linearProgramming.linearProgramming("maximize", numpy.ones(1), numpy.ones((1,1)), numpy.ones(1), numpy.ones(1), numpy.ones(1), "S", "C", None)

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/22/2013
"""

from __future__ import division

import cvxopt
import cvxopt.solvers
import numpy

def linearProgramming(maximizeFlag, f, A, b, lb, ub, constraintTypes, variableTypes, options):
	# Input validation
	if type(A) != numpy.ndarray or A.size < 1:
		raise Exception, "Matrix A must be a non-empty numpy array."
	if type(b) != numpy.ndarray or b.size != A.shape[0]:
		raise Exception, "Vector b must be a numpy array with same number of rows as A."
	if lb is None:
		lb = - numpy.ones(A.shape[1]) * numpy.Inf
	elif type(lb) != numpy.ndarray or lb.size != A.shape[1]:
		raise Exception, "Vector lb must be a numpy array with dimension equal to the number of columns of A."
	if ub is None:
		ub = numpy.ones(A.shape[1]) * numpy.Inf
	elif type(ub) != numpy.ndarray or ub.size != A.shape[1]:
		raise Exception, "Vector ub must be a numpy array with dimension equal to the number of columns of A."
	if not numpy.all(numpy.array(constraintTypes) == "S"):
		raise Exception, "At least one of your constraintTypes is not currently supported."
	if not numpy.all(numpy.array(variableTypes) == "C"):
		raise Exception, "At least one of your variableTypes is not currently supported."
	if options != None or (options != None and options["glpk"] != None):
		raise Exception, "Passing GLPK options is not currently supported."

	if maximizeFlag.lower() == "maximize":
		f = f * (-1.0)			# DO *****NOT***** DO "*=" AS THAT WILL PERFORM MODIFICATION IN PLACE

	m, n = A.shape
	G = numpy.concatenate([numpy.identity(n), -1 * numpy.identity(n)], axis = 0)
	h = numpy.concatenate([ub, -1.0 * lb], axis = 0)

	A_m = cvxopt.matrix(A)
	f_m = cvxopt.matrix(f)
	G_m = cvxopt.matrix(G)
	h_m = cvxopt.matrix(h)
	b_m = cvxopt.matrix(b)

	# Turn off screen output (http://abel.ee.ucla.edu/cvxopt/userguide/coneprog.html?highlight=solvers.lp#s-external)
	cvxopt.solvers.options["LPX_K_MSGLEV"] = 0

	sol = cvxopt.solvers.lp(f_m, G_m, h_m, A = A_m, b = b_m, solver = "glpk")

	return numpy.array(sol["x"]), sol

