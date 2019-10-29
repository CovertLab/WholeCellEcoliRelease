import time

from fireworks import FireTaskBase, explicit_serialize
from models.ecoli.sim.simulation import EcoliDaughterSimulation
from wholecell.sim.simulation import DEFAULT_SIMULATION_KWARGS

@explicit_serialize
class SimulationDaughterTask(FireTaskBase):

	_fw_name = "SimulationDaughterTask"
	required_params = ["input_sim_data", "output_directory", "inherited_state_path"]
	optional_params = ["seed",
					   "timeline",
					   "length_sec",
					   "timestep_safety_frac",
					   "timestep_max",
					   "timestep_update_freq",
					   "log_to_shell",
					   "log_to_disk_every",
					   "mass_distribution",
					   "growth_rate_noise",
					   "d_period_division",
					   "translation_supply",
					   "variable_elongation_transcription",
					   "variable_elongation_translation"]

	def run_task(self, fw_spec):

		print "%s: Running simulation" % time.ctime()

		options = {}
		
		options["simDataLocation"] = self["input_sim_data"]
		options["outputDir"] = self["output_directory"]
		options["logToDisk"] = True
		options["overwriteExistingFiles"] = False
		options["inheritedStatePath"] = self["inherited_state_path"]

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
		options["variable_elongation_transcription"] = self.get("variable_elongation_transcription", DEFAULT_SIMULATION_KWARGS["variable_elongation_transcription"])
		options["variable_elongation_translation"] = self.get("variable_elongation_translation", DEFAULT_SIMULATION_KWARGS["variable_elongation_translation"])

		sim = EcoliDaughterSimulation(**options)

		sim.run()
