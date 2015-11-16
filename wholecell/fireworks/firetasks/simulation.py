import time

from fireworks import FireTaskBase, explicit_serialize
from models.ecoli.sim.simulation import EcoliSimulation
from wholecell.sim.simulation import DEFAULT_SIMULATION_KWARGS

@explicit_serialize
class SimulationTask(FireTaskBase):

	_fw_name = "SimulationTask"
	required_params = ["input_kb", "output_directory"]
	optional_params = ["seed", "length_sec", "log_to_shell", "log_to_disk_every"]

	def run_task(self, fw_spec):

		print "%s: Running simulation" % time.ctime()

		options = {}
		
		options["kbLocation"] = self["input_kb"]
		options["outputDir"] = self["output_directory"]
		options["logToDisk"] = True
		options["overwriteExistingFiles"] = False

		options["seed"] = int(self.get("seed", DEFAULT_SIMULATION_KWARGS["seed"]))
		options["lengthSec"] = self.get("length_sec", DEFAULT_SIMULATION_KWARGS["lengthSec"])
		options["logToShell"] = self.get("log_to_shell", DEFAULT_SIMULATION_KWARGS["logToShell"])
		options["logToDiskEvery"] = self.get("log_to_disk_every", DEFAULT_SIMULATION_KWARGS["logToDiskEvery"])

		sim = EcoliSimulation(**options)

		sim.run()