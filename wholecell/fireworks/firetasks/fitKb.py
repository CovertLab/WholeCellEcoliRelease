import cPickle
import time

from fireworks import FireTaskBase, explicit_serialize
from reconstruction.ecoli.fitkb1 import fitKb_1
from reconstruction.ecoli.fitkb2 import fitKb_2

@explicit_serialize
class FitKbTask(FireTaskBase):

	_fw_name = "FitKbTask"
	required_params = ["fit_level", "input_kb", "output_kb"]
	optional_params = ["sim_out_dir"]

	def run_task(self, fw_spec):

		print "%s: Fitting knowledgebase (Level %d)" % (time.ctime(), self["fit_level"])

		if self["fit_level"] == 1:
			kb = cPickle.load(open(self["input_kb"], "rb"))
			fitKb_1(kb)
			cPickle.dump(
				kb,
				open(self["output_kb"], "wb"),
				protocol = cPickle.HIGHEST_PROTOCOL
				)

		if self["fit_level"] == 2:
			kb = cPickle.load(open(self["input_kb"], "rb"))
			fitKb_2(kb, self["sim_out_dir"])
			cPickle.dump(
				kb,
				open(self["output_kb"], "wb"),
				protocol = cPickle.HIGHEST_PROTOCOL
				)