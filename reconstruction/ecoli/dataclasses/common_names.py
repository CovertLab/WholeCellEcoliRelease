"""
Simulation Data common names.

This function creates dictionaries of common names of simulation elements in
sim_data, with {Object ID: (Primary Name, Synonyms)}.

The tsv files in 'reconstruction/ecoli/flat/common_names', which are used to
build the common_names dictionary, were taken from ecocyc and reformatted.
Ecocyc provided the columns 'Object ID' and 'Synonyms'; entries in
'Primary Name' were mostly taken from the first entry in synonyms, but in some
cases they were selected by hand.

@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import absolute_import, division, print_function

class CommonNames(object):
	""" GetNames """

	def __init__(self, raw_data, sim_data):
		self._getNamesDict(raw_data, sim_data)

	def _getNamesDict(self, raw_data, sim_data):
		# create dictionaries of common names
		names = {}

		for label in dir(raw_data.common_names):
			if label.startswith("__"):
				continue

			label_names = {}
			file_ = getattr(raw_data.common_names, label)
			for row in file_:
				label_names[row['Object ID']] = (
					row['Primary Name'],
					row['Synonyms'][1:-1].split(', '))

			names[label] = label_names

		self.__dict__.update(names)
