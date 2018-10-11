"""
Simulation Data common names.

This function creates a dictionary of common names in sim_data, with {Object ID: (Display Name, Synonyms)}.
build_network uses this dict to assign common names to the nodes in the Causality Network visualization tool.

The tsv files in 'reconstruction/ecoli/flat/common_names', which are used to build the common_names dictionary were
taken from ecocyc and reformatted. Ecocyc provided the columns 'Object ID' and 'Synonyms'; entries in 'Display Name'
were mostly taken from the first entry in synonyms, but in some cases they were selected by hand.

@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import absolute_import, division, print_function

class GetNames(object):
	""" GetNames """

	def __init__(self, raw_data, sim_data):
		self._getNamesDict(raw_data, sim_data)

	def _getNamesDict(self, raw_data, sim_data):

		# create a dictionary of all nutrient time series
		self.common_names = {}
		for label in dir(raw_data.common_names):
			if label.startswith("__"):
				continue

			file = getattr(raw_data.common_names, label)
			for row in file:
				if row['Display Name']:
					self.common_names[row['Object ID']] = (
						row['Display Name'],
						row['Synonyms']
						)
