import cPickle
import time

from fireworks import FireTaskBase, explicit_serialize
from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli

@explicit_serialize
class InitRawValidationDataTask(FireTaskBase):

	_fw_name = "InitRawValidationDataTask"
	required_params = ["output"]

	def run_task(self, fw_spec):
		print "%s: Instantiating validation_data_raw" % (time.ctime())

		validation_data_raw = ValidationDataRawEcoli()

		print "%s: Saving validation_data_raw" % (time.ctime())

		cPickle.dump(
			validation_data_raw,
			open(self["output"], "wb"),
			protocol = cPickle.HIGHEST_PROTOCOL
			)
