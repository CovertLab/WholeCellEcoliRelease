from __future__ import absolute_import, division, print_function

import cPickle
import time

from models.ecoli.sim.variants import nameToFunctionMapping


def apply_variant(sim_data_file, variant_type, variant_index):
	"""Load the sim_data_file and apply the named variant type and variant index.
	Returns (info_dict, modified_sim_data).
	Example info_dict = {'shortName': "wildtype", 'desc': "Wildtype simulation"}
	"""

	if variant_type not in nameToFunctionMapping:
		raise Exception("%s is not a valid variant function!" % variant_type)

	print("{}: Creating sim_data for Variant: {}, Index: {}".format(
		time.ctime(), variant_type, variant_index))

	with open(sim_data_file, "rb") as f:
		sim_data = cPickle.load(f)

	info, sim_data = nameToFunctionMapping[variant_type](sim_data, variant_index)

	print("Variant short name:", info["shortName"])

	return info, sim_data
