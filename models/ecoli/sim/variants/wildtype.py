"""
Wildtype variant (no changes)
Default variant for simulations

Modifies:
	nothing

Expected variant indices:
	0: wildtype
"""

from __future__ import absolute_import, division, print_function


CONTROL_OUTPUT = dict(
	shortName = "wildtype",
	desc = "Wildtype simulation"
	)


def wildtype(sim_data, index):
	return CONTROL_OUTPUT, sim_data
