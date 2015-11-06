import time

from fireworks import FireTaskBase, explicit_serialize
from models.ecoli.sim.simulation import EcoliDaughterSimulation
from wholecell.sim.simulation import DEFAULT_SIMULATION_KWARGS

@explicit_serialize
class SimulationDaughterTask(FireTaskBase):

	_fw_name = "SimulationDaughterTask"
	required_params = ["input_kb", "output_directory", "inherited_state_path"]
	optional_params = ["seed", "length_sec", "log_to_shell", "log_to_disk_every"]

	def run_task(self, fw_spec):

		print "%s: Running simulation" % time.ctime()

		options = {}
		
		options["kbLocation"] = self["input_kb"]
		options["outputDir"] = self["output_directory"]
		options["logToDisk"] = True
		options["overwriteExistingFiles"] = False
		options["inheritedStatePath"] = self["inherited_state_path"]

		options["seed"] = int(self.get("seed", DEFAULT_SIMULATION_KWARGS["seed"]))
		options["lengthSec"] = self.get("length_sec", DEFAULT_SIMULATION_KWARGS["lengthSec"])
		options["logToShell"] = self.get("log_to_shell", DEFAULT_SIMULATION_KWARGS["logToShell"])
		options["logToDiskEvery"] = self.get("log_to_disk_every", DEFAULT_SIMULATION_KWARGS["logToDiskEvery"])

		sim = EcoliDaughterSimulation(**options)

		sim.run()