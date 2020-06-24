from __future__ import absolute_import, division, print_function

from six.moves import cPickle
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

		with open(self["validation_data_input"], "rb") as data:
			raw_validation_data = cPickle.load(data)
		with open(self["knowledge_base_raw"], "rb") as raw:
			knowledge_base_raw = cPickle.load(raw)
		validation_data = ValidationDataEcoli()
		validation_data.initialize(raw_validation_data, knowledge_base_raw)

		with open(self["output_data"], "wb") as fh:
			cPickle.dump(validation_data, fh, protocol=cPickle.HIGHEST_PROTOCOL)
