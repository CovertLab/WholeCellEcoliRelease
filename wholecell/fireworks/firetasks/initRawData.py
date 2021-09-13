from __future__ import absolute_import, division, print_function

from six.moves import cPickle
import time

from fireworks import FiretaskBase, explicit_serialize
from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli
from wholecell.utils.constants import DEFAULT_OPERON_OPTION


@explicit_serialize
class InitRawDataTask(FiretaskBase):

	_fw_name = "InitRawDataTask"
	required_params = ["output"]
	optional_params = ['operons']

	def run_task(self, fw_spec):
		print("%s: Instantiating raw_data" % (time.ctime(),))

		operons_on = (self.get('operons') or DEFAULT_OPERON_OPTION) == 'on'
		raw_data = KnowledgeBaseEcoli(
			operons_on=operons_on)

		print("%s: Saving raw_data" % (time.ctime(),))

		with open(self["output"], "wb") as f:
			cPickle.dump(raw_data, f, protocol = cPickle.HIGHEST_PROTOCOL)
