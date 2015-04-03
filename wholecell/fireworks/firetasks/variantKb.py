import cPickle
import time
import os

from fireworks import FireTaskBase, explicit_serialize
from models.ecoli.sim.variants import nameToFunctionMapping

@explicit_serialize
class VariantKbTask(FireTaskBase):

	_fw_name = "VariantKbTask"
	required_params = [
		"variant_function", "variant_index",
		"input_kb", "output_kb",
		"variant_metadata_directory",
		]

	def run_task(self, fw_spec):

		if self["variant_function"] not in nameToFunctionMapping:
			raise Exception, "%s is not a valid variant function!" % self["variant_function"]

		print "%s: Creating variant kb (Variant: %s Index: %d)" % (time.ctime(), self["variant_function"], self["variant_index"])

		kb = cPickle.load(open(self["input_kb"], "rb"))

		info = nameToFunctionMapping[self["variant_function"]](kb, self["variant_index"])

		print info["shortName"]

		cPickle.dump(
			kb,
			open(self["output_kb"], "wb"),
			protocol = cPickle.HIGHEST_PROTOCOL
			)

		h = open(os.path.join(self["variant_metadata_directory"], "short_name"), "w")
		h.write("%s\n" % info["shortName"])
		h.close()

		h = open(os.path.join(self["variant_metadata_directory"], "description"), "w")
		h.write("%s\n" % info["desc"])
		h.close()