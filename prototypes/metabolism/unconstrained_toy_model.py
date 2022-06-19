from __future__ import absolute_import, division, print_function

import numpy as np
import cvxpy

metabolic_reactions = np.array([
	# R1  R2a R2b R3  R4  R5  R6  R7  R8a R8b Rres
	[-1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 ], # A 
	[ 1, -1,  1, -1,  0,  0,  0,  0,  0,  0,  0 ], # B
	[ 0,  1, -1,  0, -1, .8, -1, -1,  0,  0,  0 ], # C
	[ 0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  0 ], # D
	[ 0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0 ], # E
	[ 0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0 ], # F
	[ 0,  0,  0,  0,  1, -1,  0,  0, -1,  1,  0 ], # G
	[ 0,  0,  0,  0,  0,  0,  0,  0,  1, -1,  0 ], # H
	[-1,  2, -2,  0,  0,  0,  2,  0, -1,  1,  1 ], # Atp
	[ 0,  2, -2,  0,  0,  2,  0, -4, -2,  2, -1 ], # Nad
	[ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1 ], # O2
	[ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 ], # C1
	[ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 ], # C2
	[ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 ], # Fx
	[ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 ], # Hx
	[ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 ], # Dx
	[ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 ], # Ex
	[ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 ], # Ox
	[ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 ]  # Biomass
	])


transport_rxns = np.array([
	# Tc1 Tc2 Tf  Td  Te  Th  To2
	[ 1,  1,  0,  0,  0,  0,  0], # A
	[ 0,  0,  0,  0,  0,  0,  0], # B
	[ 0,  0,  0,  0,  0,  0,  0], # C
	[ 0,  0,  0, -1,  0,  0,  0], # D
	[ 0,  0,  0,  0, -1,  0,  0], # E
	[ 0,  0,  1,  0,  0,  0,  0], # F
	[ 0,  0,  0,  0,  0,  0,  0], # G
	[ 0,  0,  0,  0,  0,  1,  0], # H
	[ 0,  0,  0,  0,  0,  0,  0], # Atp
	[ 0,  0,  0,  0,  0,  0,  0], # Nad
	[ 0,  0,  0,  0,  0,  0,  1], # O2
	[-1,  0,  0,  0,  0,  0,  0], # C1
	[ 0, -1,  0,  0,  0,  0,  0], # C2
	[ 0,  0, -1,  0,  0,  0,  0], # Fx
	[ 0,  0,  0,  0,  0, -1,  0], # Hx
	[ 0,  0,  0,  1,  0,  0,  0], # Dx
	[ 0,  0,  0,  0,  1,  0,  0], # Ex
	[ 0,  0,  0,  0,  0,  0, -1], # Ox
	[ 0,  0,  0,  0,  0,  0,  0]  # Biomass
])


source_reactions = np.array([
	# C1  C2  Fx  Hx  Dx  Ex  Ox  Growth
	[ 0,  0,  0,  0,  0,  0,  0,  0],  # A
	[ 0,  0,  0,  0,  0,  0,  0,  0],  # B
	[ 0,  0,  0,  0,  0,  0,  0,  0],  # C
	[ 0,  0,  0,  0,  0,  0,  0,  0],  # D
	[ 0,  0,  0,  0,  0,  0,  0,  0],  # E
	[ 0,  0,  0,  0,  0,  0,  0,  0],  # F
	[ 0,  0,  0,  0,  0,  0,  0,  0],  # G
	[ 0,  0,  0,  0,  0,  0,  0,  0],  # H
	[ 0,  0,  0,  0,  0,  0,  0,  0],  # Atp
	[ 0,  0,  0,  0,  0,  0,  0,  0],  # Nad
	[ 0,  0,  0,  0,  0,  0,  0,  0],  # O2
	[ 1,  0,  0,  0,  0,  0,  0,  0],  # C1
	[ 0,  1,  0,  0,  0,  0,  0,  0],  # C2
	[ 0,  0,  1,  0,  0,  0,  0,  0],  # Fx
	[ 0,  0,  0,  1,  0,  0,  0,  0],  # Hx
	[ 0,  0,  0,  0, -1,  0,  0,  0],  # Dx
	[ 0,  0,  0,  0,  0, -1,  0,  0],  # Ex
	[ 0,  0,  0,  0,  0,  0,  1,  0],  # Ox
	[ 0,  0,  0,  0,  0,  0,  0, -1]   # Biomass
	])

biomass_reaction = np.array([
	# Biomass reaction
	[  0],  # A
	[  0],  # B
	[ -1],  # C
	[  0],  # D
	[  0],  # E
	[ -1],  # F
	[  0],  # G
	[ -1],  # H
	[-10],  # Atp
	[  0],  # Nad
	[  0],  # O2
	[  0],  # C1
	[  0],  # C2
	[  0],  # Fx
	[  0],  # Hx
	[  0],  # Dx
	[  0],  # Ex
	[  0],  # Ox
	[  1]   # Biomass
	])

max_fluxes = np.array([
	[np.inf], # R1
	[np.inf], # R2a
	[np.inf], # R2b
	[np.inf], # R3
	[np.inf], # R4
	[np.inf], # R5
	[np.inf], # R6
	[np.inf], # R7
	[np.inf], # R8a
	[np.inf], # R8b
	[np.inf], # Rres

	[10.5], # C1
	[10.5], # C2
	[5.0],  # F
	[12.0], # D
	[12.0], # E
	[5.0],  # H
	[15.0], # O2

	[np.inf], # C1
	[np.inf], # C2
	[np.inf], # Fx
	[np.inf], # Hx
	[np.inf], # Dx
	[np.inf], # Ex
	[np.inf], # Ox
	[np.inf], # Growth

	[np.inf], # Biomass
	])

min_fluxes = np.array([
	[0.], # R1
	[0.], # R2a
	[0.], # R2b
	[0.], # R3
	[0.], # R4
	[0.], # R5
	[0.], # R6
	[0.], # R7
	[0.], # R8a
	[0.], # R8b
	[0.], # Rres

	[0.], # C1
	[0.], # C2
	[0.], # F
	[0.], # D
	[0.], # E
	[0.], # H
	[0.], # O2

	[0.], # C1
	[0.], # C2
	[0.], # Fx
	[0.], # Hx
	[0.], # Dx
	[0.], # Ex
	[0.], # Ox
	[0.], # Growth

	[0.], # Biomass
	])


S_matrix = np.hstack((metabolic_reactions,  transport_rxns, source_reactions,  biomass_reaction))

num_metabolites, num_reactions = S_matrix.shape

# Construct the problem
fluxes = cvxpy.Variable(num_reactions)

# One-hot vector indicating location of biomass reaction
c = np.zeros(num_reactions)
c[-1] = 1

# Maximize biomass reaction
objective = cvxpy.Maximize(c*fluxes)

constraints = [
	S_matrix*fluxes == 0,
	min_fluxes <= fluxes,
	fluxes <= max_fluxes,
	]

prob = cvxpy.Problem(objective, constraints)

# The optimal objective is returned by prob.solve().
result = prob.solve(solver=cvxpy.GUROBI)

print(fluxes.value)
