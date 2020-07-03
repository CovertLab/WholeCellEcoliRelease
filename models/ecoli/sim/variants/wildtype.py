"""
Wildtype variant (no changes)
Default variant for simulations

Modifies:
	nothing

Expected variant indices:
	0: wildtype
	1: wildtype  # no changes
	*: wildtype  # no changes

Running runSim.py with `-v wildtype 0 1` gives two variants with identical output.
This is useful for testing repeatability.
"""

from __future__ import absolute_import, division, print_function


CONTROL_OUTPUT = dict(
	shortName = "wildtype",
	desc = "Wildtype simulation"
	)


def wildtype(sim_data, index):
	_ = index
	return CONTROL_OUTPUT, sim_data
