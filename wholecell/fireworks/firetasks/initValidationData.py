from __future__ import absolute_import, division, print_function

import cPickle
import time

from fireworks import FiretaskBase, explicit_serialize
from validation.ecoli.validation_data import ValidationDataEcoli


@explicit_serialize
class InitValidationDataTask(FiretaskBase):

	_fw_name = "InitValidationDataTask"
	required_params = [
		"validation_data_input",
		"knowledge_base_raw",
		"output_data"]

	def run_task(self, fw_spec):
		print("{}: Initializing Validation Data".format(time.ctime()))

		raw_validation_data = cPickle.load(open(self["validation_data_input"], "rb"))
		knowledge_base_raw = cPickle.load(open(self["knowledge_base_raw"], "rb"))
		validation_data = ValidationDataEcoli()
		validation_data.initialize(raw_validation_data, knowledge_base_raw)
		cPickle.dump(
			validation_data,
			open(self["output_data"], "wb"),
			protocol = cPickle.HIGHEST_PROTOCOL
			)
