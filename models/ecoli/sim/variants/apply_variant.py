from __future__ import absolute_import, division, print_function

from six.moves import cPickle
import time

from models.ecoli.sim.variants import nameToFunctionMapping


def apply_variant(sim_data_file, variant_type, variant_index):
	"""Load the sim_data_file and apply the named variant type and variant index.
	Returns (info_dict, modified_sim_data).
	Example info_dict = {'shortName': "wildtype", 'desc': "Wildtype simulation"}
	"""

	if variant_type not in nameToFunctionMapping:
		raise Exception("%s is not a valid variant function!" % variant_type)

	with open(sim_data_file, "rb") as f:
		sim_data = cPickle.load(f)

	operon_msg = (f", Operons: {'on' if sim_data.operons_on else 'off'}"
				  if hasattr(sim_data, 'operons_on') else '')
	print(f"{time.ctime()}: Creating sim_data for Variant: {variant_type},"
		  f" Index: {variant_index}{operon_msg}")

	info, sim_data = nameToFunctionMapping[variant_type](sim_data, variant_index)

	print("Variant short name:", info["shortName"])

	return info, sim_data
