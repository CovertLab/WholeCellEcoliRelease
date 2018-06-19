import cPickle
import time

from fireworks import FireTaskBase, explicit_serialize
from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli

@explicit_serialize
class InitRawDataTask(FireTaskBase):

	_fw_name = "InitRawDataTask"
	required_params = ["output"]

	def run_task(self, fw_spec):
		print "%s: Instantiating raw_data" % (time.ctime(),)

		raw_data = KnowledgeBaseEcoli()

		print "%s: Saving raw_data" % (time.ctime(),)

		with open(self["output"], "wb") as f:
			cPickle.dump(raw_data, f, protocol = cPickle.HIGHEST_PROTOCOL)
