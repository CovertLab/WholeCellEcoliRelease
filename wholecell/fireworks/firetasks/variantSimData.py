from __future__ import absolute_import, division, print_function

import cPickle
import os
import sys

from fireworks import FireTaskBase, explicit_serialize
from models.ecoli.sim.variants import variant


@explicit_serialize
class VariantSimDataTask(FireTaskBase):

	_fw_name = "VariantSimDataTask"
	required_params = [
		"variant_function", "variant_index",
		"input_sim_data", "output_sim_data",
		"variant_metadata_directory",
		]

	def run_task(self, fw_spec):
		info, sim_data = variant.apply_variant(self["input_sim_data"], self["variant_function"], self["variant_index"])

		sys.setrecursionlimit(4000)

		with open(self["output_sim_data"], "wb") as f:
			cPickle.dump(sim_data, f, protocol = cPickle.HIGHEST_PROTOCOL)

		with open(os.path.join(self["variant_metadata_directory"], "short_name"), "w") as h:
			h.write("%s\n" % info["shortName"])

		with open(os.path.join(self["variant_metadata_directory"], "description"), "w") as h:
			h.write("%s\n" % info["desc"])
