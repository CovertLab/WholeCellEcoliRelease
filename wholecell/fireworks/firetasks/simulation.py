import time

from fireworks import FireTaskBase, explicit_serialize
from models.ecoli.sim.simulation import EcoliSimulation
from wholecell.sim.simulation import DEFAULT_SIMULATION_KWARGS

@explicit_serialize
class SimulationTask(FireTaskBase):

	_fw_name = "SimulationTask"
	required_params = ["input_sim_data", "output_directory"]
	optional_params = ["seed", "length_sec", "timestep_safety_frac", "timestep_max", "timestep_update_freq", "log_to_shell", "log_to_disk_every", "mass_distribution", "growth_rate_noise", "d_period_division", "translation_supply"]

	def run_task(self, fw_spec):

		print "%s: Running simulation" % time.ctime()

		options = {}
		
		options["simDataLocation"] = self["input_sim_data"]
		options["outputDir"] = self["output_directory"]
		options["logToDisk"] = True
		options["overwriteExistingFiles"] = False

		options["seed"] = int(self.get("seed", DEFAULT_SIMULATION_KWARGS["seed"]))
		options["lengthSec"] = self.get("length_sec", DEFAULT_SIMULATION_KWARGS["lengthSec"])
		options["timeStepSafetyFraction"] = self.get("timestep_safety_frac", DEFAULT_SIMULATION_KWARGS["timeStepSafetyFraction"])
		options["maxTimeStep"] = self.get("timestep_max", DEFAULT_SIMULATION_KWARGS["maxTimeStep"])
		options["updateTimeStepFreq"] = self.get("timestep_update_freq", DEFAULT_SIMULATION_KWARGS["updateTimeStepFreq"])
		options["logToShell"] = self.get("log_to_shell", DEFAULT_SIMULATION_KWARGS["logToShell"])
		options["logToDiskEvery"] = self.get("log_to_disk_every", DEFAULT_SIMULATION_KWARGS["logToDiskEvery"])
		options["massDistribution"] = self.get("mass_distribution", DEFAULT_SIMULATION_KWARGS["massDistribution"])
		options["growthRateNoise"] = self.get("growth_rate_noise", DEFAULT_SIMULATION_KWARGS["growthRateNoise"])
		options["dPeriodDivision"] = self.get("d_period_division", DEFAULT_SIMULATION_KWARGS["dPeriodDivision"])
		options["translationSupply"] = self.get("translation_supply", DEFAULT_SIMULATION_KWARGS["translationSupply"])

		sim = EcoliSimulation(**options)

		sim.run()