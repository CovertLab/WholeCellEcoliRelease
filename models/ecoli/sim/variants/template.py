"""
UPDATE: describe what variant should do and any analysis plots to use with it
UPDATE: change Modifies and Expected variant indices sections below (examples shown can be removed)

Modifies:
	sim_data.condition
	sim_data.external_state.current_timeline_id

Expected variant indices (dependent on sorted order of sim_data.conditionActiveTfs):
	0: control
	1: anaerobic
	2: with amino acids
"""

from __future__ import absolute_import, division, print_function


# UPDATE: give a descriptive function name that matches the file name
def template(sim_data, index):
	# UPDATE: modify sim_data attributes based on the variant index
	sim_data.something = index

	# UPDATE: set strings to give a description of the variant/what has changed
	return dict(
		shortName = '{}'.format(index),
		desc = 'something is set at {}.'.format(index)
		), sim_data
