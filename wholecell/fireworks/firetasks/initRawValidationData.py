from __future__ import absolute_import, division, print_function

from six.moves import cPickle
import time

from fireworks import FiretaskBase, explicit_serialize
from validation.ecoli.validation_data_raw import ValidationDataRawEcoli


@explicit_serialize
class InitRawValidationDataTask(FiretaskBase):

	_fw_name = "InitRawValidationDataTask"
	required_params = ["output"]

	def run_task(self, fw_spec):
		print("%s: Instantiating validation_data_raw" % (time.ctime()))

		validation_data_raw = ValidationDataRawEcoli()

		print("%s: Saving validation_data_raw" % (time.ctime()))

		cPickle.dump(
			validation_data_raw,
			open(self["output"], "wb"),
			protocol = cPickle.HIGHEST_PROTOCOL
			)
