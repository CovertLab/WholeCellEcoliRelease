import numpy as np
import cvxpy

# metabolic_reactions = np.array([
# 	# A   B   C   D   E   F   G   H  Atp Nad  O2  C1  C2  Fx  Hx  Dx  Ex  Ox  Biomass
# 	[-1,  1,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0], # R1
# 	[ 0, -1,  1,  0,  0,  0,  0,  0,  2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0], # R2
# 	[ 0, -1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0], # R3
# 	[ 0,  0, -1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0], # R4
# 	[ 0,  0, .8,  0,  0,  0, -1,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0], # R5
# 	[ 0,  0, -1,  3,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0], # R6
# 	[ 0,  0, -1,  0,  3,  0,  0,  0,  0, -4,  0,  0,  0,  0,  0,  0,  0,  0,  0], # R7
# 	[ 0,  0,  0,  0,  0,  0, -1,  1, -1, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0], # R8
# 	[ 0,  0,  0,  0,  0,  0,  0,  0,  1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0], # Rres
# 	])

# transport_rxns = np.array([
# 	# A   B   C   D   E   F   G   H  Atp Nad  O2  C1  C2  Fx  Hx  Dx  Ex  Ox  Biomass
# 	[ 1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0], # Tc1
# 	[ 0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0], # Tf
# 	[ 0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0], # Td
# 	[ 0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0], # Te
# 	[ 0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0], # Th
# 	[ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  1,  0], # To2
# 	])

# source_reactions = np.array([
# 	# A   B   C   D   E   F   G   H  Atp Nad  O2  C1  C2  Fx  Hx  Dx  Ex  Ox  Biomass
# 	[ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0], # C1
# 	[ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0], # C2
# 	[ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0], # Fx
# 	[ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0], # Hx
# 	[ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0], # Dx
# 	[ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0], # Ex
# 	[ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0], # Ox
# 	])

# biomass_reaction = np.array([
# 	# A   B   C   D   E   F   G   H  Atp Nad  O2  C1  C2  Fx  Hx  Dx  Ex  Ox  Biomass
# 	[ 0,  0, -1,  0,  0, -1,  0, -1,-10,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1], 
# 	])

metabolic_reactions = np.array([
	# R1  R2  R3  R4  R5  R6  R7  R8  Rres
	[-1,  0,  0,  0,  0,  0,  0,  0,  0 ], # A 
	[ 1, -1, -1,  0,  0,  0,  0,  0,  0 ], # B
	[ 0,  1,  0, -1, .8, -1, -1,  0,  0 ], # C
	[ 0,  0,  0,  0,  0,  3,  0,  0,  0 ], # D
	[ 0,  0,  0,  0,  0,  0,  3,  0,  0 ], # E
	[ 0,  0,  1,  0,  0,  0,  0,  0,  0 ], # F
	[ 0,  0,  0,  1, -1,  0,  0, -1,  0 ], # G
	[ 0,  0,  0,  0,  0,  0,  0,  1,  0 ], # H
	[-1,  2,  0,  0,  0,  2,  0, -1,  1 ], # Atp
	[ 0,  2,  0,  0,  2,  0, -4, -2, -1 ], # Nad
	[ 0,  0,  0,  0,  0,  0,  0,  0, -1 ], # O2
	[ 0,  0,  0,  0,  0,  0,  0,  0,  0 ], # C1
	[ 0,  0,  0,  0,  0,  0,  0,  0,  0 ], # C2
	[ 0,  0,  0,  0,  0,  0,  0,  0,  0 ], # Fx
	[ 0,  0,  0,  0,  0,  0,  0,  0,  0 ], # Hx
	[ 0,  0,  0,  0,  0,  0,  0,  0,  0 ], # Dx
	[ 0,  0,  0,  0,  0,  0,  0,  0,  0 ], # Ex
	[ 0,  0,  0,  0,  0,  0,  0,  0,  0 ], # Ox
	# [ 0,  0,  0,  0,  0,  0,  0,  0,  0 ]  # Biomass
	])

# transport_rxns = np.array([
# 	# Tc1 Tf  Td  Te  Th  To2
# 	[ 1,  0,  0,  0,  0,  0], # A
# 	[ 0,  0,  0,  0,  0,  0], # B
# 	[ 0,  0,  0,  0,  0,  0], # C
# 	[ 0,  0, -1,  0,  0,  0], # D
# 	[ 0,  0,  0, -1,  0,  0], # E
# 	[ 0,  1,  0,  0,  0,  0], # F
# 	[ 0,  0,  0,  0,  0,  0], # G
# 	[ 0,  0,  0,  0,  1,  0], # H
# 	[ 0,  0,  0,  0,  0,  0], # Atp
# 	[ 0,  0,  0,  0,  0,  0], # Nad
# 	[ 0,  0,  0,  0,  0,  1], # O2
# 	[-1,  0,  0,  0,  0,  0], # C1
# 	[ 0,  0,  0,  0,  0,  0], # C2
# 	[ 0, -1,  0,  0,  0,  0], # Fx
# 	[ 0,  0,  0,  0, -1,  0], # Hx
# 	[ 0,  0,  1,  0,  0,  0], # Dx
# 	[ 0,  0,  0,  1,  0,  0], # Ex
# 	[ 0,  0,  0,  0,  0, -1], # Ox
# 	[ 0,  0,  0,  0,  0,  0]  # Biomass
# ])

transport_rxns = np.array([
	# Tc1 Tf  Td  Te  Th  To2
	[ 1,  0,  0,  0,  0,  0], # A
	[ 0,  0,  0,  0,  0,  0], # B
	[ 0,  0,  0,  0,  0,  0], # C
	[ 0,  0, -1,  0,  0,  0], # D
	[ 0,  0,  0, -1,  0,  0], # E
	[ 0,  1,  0,  0,  0,  0], # F
	[ 0,  0,  0,  0,  0,  0], # G
	[ 0,  0,  0,  0,  1,  0], # H
	[ 0,  0,  0,  0,  0,  0], # Atp
	[ 0,  0,  0,  0,  0,  0], # Nad
	[ 0,  0,  0,  0,  0,  1], # O2
	[ 0,  0,  0,  0,  0,  0], # C1
	[ 0,  0,  0,  0,  0,  0], # C2
	[ 0,  0,  0,  0,  0,  0], # Fx
	[ 0,  0,  0,  0,  0,  0], # Hx
	[ 0,  0,  0,  0,  0,  0], # Dx
	[ 0,  0,  0,  0,  0,  0], # Ex
	[ 0,  0,  0,  0,  0,  0], # Ox
	# [ 0,  0,  0,  0,  0,  0]  # Biomass
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
	[ 0,  0,  0,  0,  1,  0,  0,  0],  # Dx
	[ 0,  0,  0,  0,  0,  1,  0,  0],  # Ex
	[ 0,  0,  0,  0,  0,  0,  1,  0],  # Ox
	# [ 0,  0,  0,  0,  0,  0,  0, -1]   # Biomass
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
	# [  0]   # Biomass
	])

max_fluxes = np.array([
	[np.inf], # R1
	[np.inf], # R2
	[np.inf], # R3
	[np.inf], # R4
	[np.inf], # R5
	[np.inf], # R6
	[np.inf], # R7
	[np.inf], # R8
	[np.inf], # Rres

	[10.5], # C1
	[5.0],  # F
	[12.0], # D
	[12.0], # E
	[5.0],  # H
	[15.0], # O2

	# [np.inf], # C1
	# [np.inf], # C2
	# [np.inf], # Fx
	# [np.inf], # Hx
	# [np.inf], # Dx
	# [np.inf], # Ex
	# [np.inf], # Ox
	# [np.inf], # Growth

	[np.inf], # Biomass
	])


min_fluxes = np.array([
	[-np.inf], # R1
	[-np.inf], # R2
	[-np.inf], # R3
	[-np.inf], # R4
	[-np.inf], # R5
	[-np.inf], # R6
	[-np.inf], # R7
	[-np.inf], # R8
	[-np.inf], # Rres

	[0.], # C1
	[0.], # F
	[0.], # D
	[0.], # E
	[0.], # H
	[0.], # O2

	# [-np.inf], # C1
	# [-np.inf], # C2
	# [-np.inf], # Fx
	# [-np.inf], # Hx
	# [-np.inf], # Dx
	# [-np.inf], # Ex
	# [-np.inf], # Ox
	# [-np.inf], # Growth

	[-np.inf], # Biomass
	])


# S_matrix = np.hstack((metabolic_reactions,  transport_rxns, source_reactions,  biomass_reaction))
S_matrix = np.hstack((metabolic_reactions,  transport_rxns, biomass_reaction))


num_metabolites, num_reactions = S_matrix.shape

# Construct the problem.
fluxes = cvxpy.Variable(num_reactions)

# One-hot vector indicating location of biomass reaction
# c = np.zeros((1,num_reactions))
# c[:,-1] = 1

# Maximize biomass reaction
objective = cvxpy.Maximize(fluxes[num_reactions-1])

constraints = [
	S_matrix*fluxes == 0,
	min_fluxes.reshape(-1,1) <= fluxes,
	fluxes <= np.inf*np.ones_like(max_fluxes.reshape(-1,1))
	]

prob = cvxpy.Problem(objective, constraints)

# The optimal objective is returned by prob.solve().
result = prob.solve()
# The optimal value for x is stored in x.value.
print fluxes.value
# The optimal Lagrange multiplier for a constraint
# is stored in constraint.dual_value.
# print constraints[0].dual_value