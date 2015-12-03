import cPickle
import time

from fireworks import FireTaskBase, explicit_serialize
from validation.ecoli.validation_data import ValidationDataEcoli

@explicit_serialize
class InitValidationDataTask(FireTaskBase):

	_fw_name = "InitValidationDataTask"
	required_params = ["input_data", "output_data"]

	def run_task(self, fw_spec):

		print "%s: Initializing Validation Data" % (time.ctime())

		raw_data = cPickle.load(open(self["input_data"], "rb"))
		validation_data = ValidationDataEcoli()
		validation_data._initialize(raw_data)
		cPickle.dump(
			validation_data,
			open(self["output_data"], "wb"),
			protocol = cPickle.HIGHEST_PROTOCOL
			)