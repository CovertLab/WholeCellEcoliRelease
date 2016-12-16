import cPickle
import time
import os

from fireworks import FireTaskBase, explicit_serialize
from models.ecoli.sim.variants import nameToFunctionMapping

@explicit_serialize
class VariantSimDataTask(FireTaskBase):

	_fw_name = "VariantSimDataTask"
	required_params = [
		"variant_function", "variant_index",
		"input_sim_data", "output_sim_data",
		"variant_metadata_directory",
		]

	def run_task(self, fw_spec):

		if self["variant_function"] not in nameToFunctionMapping:
			raise Exception, "%s is not a valid variant function!" % self["variant_function"]

		print "%s: Creating variant sim_data (Variant: %s Index: %d)" % (time.ctime(), self["variant_function"], self["variant_index"])

		sim_data = cPickle.load(open(self["input_sim_data"], "rb"))

		info, sim_data = nameToFunctionMapping[self["variant_function"]](sim_data, self["variant_index"])

		print info["shortName"]

		import sys; sys.setrecursionlimit(4000)

		cPickle.dump(
			sim_data,
			open(self["output_sim_data"], "wb"),
			protocol = cPickle.HIGHEST_PROTOCOL
			)

		h = open(os.path.join(self["variant_metadata_directory"], "short_name"), "w")
		h.write("%s\n" % info["shortName"])
		h.close()

		h = open(os.path.join(self["variant_metadata_directory"], "description"), "w")
		h.write("%s\n" % info["desc"])
		h.close()