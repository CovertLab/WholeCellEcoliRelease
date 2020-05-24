from __future__ import absolute_import, division, print_function

import numpy as np
import cvxpy

def main():

	metabolic_reactions = np.array([
		# R1  R2
		[-1,  0], # A
		[-2,  0], # B
		[ 1, -2], # C
		[ 0,  1], # D
		[ 0,  0], # Aext
		[ 0,  0], # Bext
		[ 0,  0], # biomass
		])

	transport_rxns = np.array([
		# T1  T2
		[ 1,  0], # A
		[ 0,  1], # B
		[ 0,  0], # C
		[ 0,  0], # D
		[-1,  0], # Aext
		[ 0, -1], # Bext
		[ 0,  0], # biomass
		])

	source_reactions = np.array([
		# SA  SB  Sbiomass
		[ 0,  0,  0], # A
		[ 0,  0,  0], # B
		[ 0,  0,  0], # C
		[ 0,  0,  0], # D
		[ 1,  0,  0], # Aext
		[ 0,  1,  0], # Bext
		[ 0,  0, -1], # biomass
		])

	biomass_reaction = np.array([
		[ 0], # A
		[ 0], # B
		[-2], # C
		[-10], # D
		[ 0], # Aext
		[ 0], # Bext
		[ 1], # biomass
		])

	max_fluxes = np.array([
		[np.inf], # R1
		[np.inf], # R2
		[10.], # T1
		[5], # T2
		[np.inf], # SA
		[np.inf], # SB
		[np.inf], # Sbiomass
		[np.inf], # Biomass
		])


	min_fluxes = np.array([
		[0], # R1
		[0], # R2
		[0], # T1
		[0], # T2
		[0], # SA
		[0], # SB
		[0], # Sbiomass
		[0], # Biomass
		])

	# metabolic_reactions = np.array([
	# 	# R1
	# 	[-1], # A
	# 	[-2], # B
	# 	[ 0], # Aext
	# 	[ 0], # biomass
	# 	])

	# transport_rxns = np.array([
	# 	# T1
	# 	[-1], # A
	# 	[-2], # B
	# 	[ 0], # Aext
	# 	[ 0], # biomass
	# 	])

	# source_reactions = np.array([
	# 	# SA  Sbiomass
	# 	[ 0,  0], # A
	# 	[ 0,  0], # B
	# 	[ 1,  0], # Aext
	# 	[ 0, -1], # biomass
	# 	])

	# biomass_reaction = np.array([
	# 	[ 0], # A
	# 	[-1], # B
	# 	[ 0], # Aext
	# 	[ 1], # biomass
	# 	])


	# max_fluxes = np.array([
	# 	[np.inf], # R1
	# 	[np.inf], # T1
	# 	[np.inf], # SA
	# 	[np.inf], # Sbiomass
	# 	[np.inf], # Biomass
	# 	])


	# min_fluxes = np.array([
	# 	[-np.inf], # R1
	# 	[-np.inf], # T1
	# 	[-np.inf], # SA
	# 	[-np.inf], # Sbiomass
	# 	[-np.inf], # Biomass
	# 	])


	# S_matrix = np.hstack((metabolic_reactions,  transport_rxns, source_reactions,  biomass_reaction))
	S_matrix = np.hstack((metabolic_reactions, transport_rxns, source_reactions, biomass_reaction))

	num_metabolites, num_reactions = S_matrix.shape

	# Construct the problem.
	fluxes = cvxpy.Variable(num_reactions)

	# One-hot vector indicating location of biomass reaction
	c = np.zeros((1,num_reactions))
	c[:,-1] = 1

	# Maximize biomass reaction
	objective = cvxpy.Maximize(c*fluxes)

	constraints = [
		S_matrix*fluxes == 0,
		min_fluxes.reshape(-1,1) <= fluxes,
		fluxes <= max_fluxes.reshape(-1,1)
		]

	prob = cvxpy.Problem(objective, constraints)

	# The optimal objective is returned by prob.solve().
	result = prob.solve(solver=cvxpy.GUROBI)
	# The optimal value for x is stored in x.value.
	print(fluxes.value)
	# The optimal Lagrange multiplier for a constraint
	# is stored in constraint.dual_value.
	# print(constraints[0].dual_value)

if __name__ == "__main__":
	main()
