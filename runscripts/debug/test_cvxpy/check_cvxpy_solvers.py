'''
Test non-deterministic behavior between local and sherlock parca outputs.
Includes necessary pieces to reproduce from fitPromoterBoundProbability
function from fit_sim_data_1.py.

Results for solver options:
  GLPK (original method): differences locally and in docker vs sherlock but gives 1.0 and 0.0 when expected
  ECOS: gives comparable results and is faster (6-20x) but floating point issues for 0 and 1
  SCS: gives comparable results and is faster (1-3x) but has floating point issues for 0 and 1
  OSQP: gives different results and has floating point issues for 0 and 1 but is faster
'''

import os
import time

from cvxpy import Variable, Problem, Minimize, norm
import numpy as np


solvers = ['GLPK', 'ECOS', 'SCS', 'OSQP']
cached_dir = os.path.join(os.path.dirname(__file__), 'saved_arrays')

# Parameters used in fitPromoterBoundProbability()
PROMOTER_PDIFF_THRESHOLD = 0.1  # Minimum difference between binding probabilities of a TF in conditions where TF is active and inactive
PROMOTER_REG_COEFF = 1e-3  # Optimization weight on how much probability should stay close to original values
PROMOTER_SCALING = 10  # Multiplied to all matrices for numerical stability
PROMOTER_NORM_TYPE = 1  # Matrix 1-norm


def test_solver(problem, output, solver):
	'''
	Print output for each method along with time to solve

	Args:
		problem (cvxpy.Problem object): problem setup
		output (cvxpy.Variable object): output variable to check
		solver (str): solver method to use
	'''

	start = time.time()

	try:
		problem.solve(solver=solver)
	except Exception as e:
		print('\nSolving with {} failed: {}'.format(solver, e))
		return

	print('\n{} solved in {:.2f} s'.format(solver, time.time() - start))

	sol = np.array(output.value).reshape(-1)

	# Example indices for desired outputs, other indices can be checked
	print(sol[114])  # Might be different locally vs sherlock with GLPK
	print(sol[0])  # Should be 1.0
	print(sol[21])  # Should be 0.0


# Reproduced from fitPromoterBoundProbability()
k = np.load(os.path.join(cached_dir, 'k.npy'))

# Use optimal value of R to construct matrix H and vector Pdiff
Hi = np.load(os.path.join(cached_dir, 'Hi.npy'))
Hj = np.load(os.path.join(cached_dir, 'Hj.npy'))
Hv = np.load(os.path.join(cached_dir, 'Hv.npy'))
H = np.zeros((Hi.max() + 1, Hj.max() + 1))
H[Hi, Hj] = Hv
pInit0 = np.load(os.path.join(cached_dir, 'pInit0.npy'))
pdiff = np.load(os.path.join(cached_dir, 'pdiff.npy'))

D = np.load(os.path.join(cached_dir, 'D.npy'))
Drhs = np.load(os.path.join(cached_dir, 'Drhs.npy'))

# Optimize P such that the transcription probabilities computed from
# current values of R in matrix H are close to fit values.
P = Variable(H.shape[1])

# Objective: minimize difference between k (fit RNAP initiation
# probabilities) and H*P (computed initiation probabilities) while
# also minimizing deviation of P from the original value calculated
# from mean TF and ligand concentrations
objective_p = Minimize(
	norm(H * (PROMOTER_SCALING * P) - PROMOTER_SCALING * k, PROMOTER_NORM_TYPE)
	+ PROMOTER_REG_COEFF * norm(P - pInit0, PROMOTER_NORM_TYPE)
	)

# Constraints
# 1) 0 <= P <= 1 : All DNA-bound probabilities should be between zero
# and one.
# 2) D*P == Drhs : Values of P that correspond to alpha's and fixed TFs
# should not change.
# 3) pdiff*P >= 0.1 : There must be at least a difference of 0.1
# between binding probabilities of a TF in conditions TF__active and
# TF__inactive
constraint_p = [
	0 <= PROMOTER_SCALING*P, PROMOTER_SCALING*P <= PROMOTER_SCALING*1,
	np.diag(D)*(PROMOTER_SCALING*P) == PROMOTER_SCALING*Drhs,
	pdiff*(PROMOTER_SCALING*P) >= PROMOTER_SCALING*PROMOTER_PDIFF_THRESHOLD
	]

# Solve optimization problem
prob_p = Problem(objective_p, constraint_p)

for solver in solvers:
	test_solver(prob_p, P, solver)
