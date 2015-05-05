import cPickle
import time

from fireworks import FireTaskBase, explicit_serialize
import reconstruction.ecoli.simulation_data

@explicit_serialize
class InitKbTask(FireTaskBase):

	_fw_name = "InitKbTask"
	required_params = ["output"]

	def run_task(self, fw_spec):
		print "%s: Instantiating unfit knowledgebase" % (time.ctime())

		kb = reconstruction.ecoli.simulation_data.SimulationDataEcoli()

		print "%s: Saving unfit knowledgebase" % (time.ctime())

		cPickle.dump(
			kb,
			open(self["output"], "wb"),
			protocol = cPickle.HIGHEST_PROTOCOL
			)
