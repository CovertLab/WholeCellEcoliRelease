"""
Base implementation of any solver to be used with the parameter_search.py
runscript.  Will update parameters (raw or sim data), run the parca and run
sims to get an objective to optimize for.  Specific implementations should
subclass from BaseSolver and implement the required functions.
"""

import abc
import os
import pickle
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from models.ecoli.sim.variants.apply_variant import apply_variant
from models.ecoli.sim.parameter_search.base_parameter_search import Parameter
from wholecell.fireworks.firetasks import FitSimDataTask, InitRawDataTask, SimulationTask, SimulationDaughterTask, VariantSimDataTask
from wholecell.sim.simulation import ALTERNATE_KWARG_NAMES
from wholecell.utils import constants, data, parallelization, scriptBase, units
import wholecell.utils.filepath as fp


def run_sim(args):
	for alt, original in ALTERNATE_KWARG_NAMES.items():
		if alt not in args:
			args[alt] = args[original]

	sim_args = data.select_keys(args, scriptBase.SIM_KEYS)

	variant_directory = args['variant directory']
	variant_sim_data_directory = os.path.join(variant_directory, constants.VKB_DIR)
	variant_sim_data_modified_file = os.path.join(variant_sim_data_directory, constants.SERIALIZED_SIM_DATA_MODIFIED)

	# TODO: should this be applied before sim_data parameters are updated?
	if 'variant' in args:
		variant, index = args['variant']
		_, sim_data = apply_variant(variant_sim_data_modified_file, variant, index)
		updated_sim_data_file = os.path.join(
			os.path.dirname(variant_sim_data_modified_file),
			f'{args["index"]}-{os.path.basename(variant_sim_data_modified_file)}',
			)

		with open(updated_sim_data_file, 'wb') as f:
			pickle.dump(sim_data, f, protocol=pickle.HIGHEST_PROTOCOL)

		variant_sim_data_modified_file = updated_sim_data_file

	for j in range(args['seed'], args['seed'] + args['init_sims']):  # init sim seeds, TODO: allow for multiple seeds with index
		seed_directory = fp.makedirs(variant_directory, "%06d" % args['index'])

		for k in range(args['generations']):  # generation number k
			gen_directory = fp.makedirs(seed_directory, "generation_%06d" % k)

			# l is the daughter number among all of this generation's cells,
			# which is 0 for single-daughters but would span range(2**k) if
			# each parent had 2 daughters.
			l = 0
			cell_directory = fp.makedirs(gen_directory, "%06d" % l)
			cell_sim_out_directory = fp.makedirs(cell_directory, "simOut")

			options = dict(sim_args,
				input_sim_data=variant_sim_data_modified_file,
				output_directory=cell_sim_out_directory,
				)

			if k == 0:
				task = SimulationTask(seed=j, **options)
			else:
				parent_gen_directory = os.path.join(seed_directory, "generation_%06d" % (k - 1))
				parent_cell_directory = os.path.join(parent_gen_directory, "%06d" % (l // 2))
				parent_cell_sim_out_directory = os.path.join(parent_cell_directory, "simOut")
				daughter_state_path = os.path.join(
					parent_cell_sim_out_directory,
					constants.SERIALIZED_INHERITED_STATE % (l % 2 + 1))
				task = SimulationDaughterTask(
					seed=(j + 1) * ((2 ** k - 1) + l),
					inherited_state_path=daughter_state_path,
					**options
					)
			task.run_task({})

	return variant_directory

def run_parca(args, raw_data_file, sim_data_file, metrics_data_file, cpus):
	task = FitSimDataTask(
		input_data=raw_data_file,
		output_data=sim_data_file,
		output_metrics_data=metrics_data_file,
		cached=False,
		load_intermediate=None,
		save_intermediates=False,
		intermediates_directory=os.path.dirname(sim_data_file),
		cpus=cpus,
		debug=args.get('debug', False),
		variable_elongation_transcription=args.get('variable_elongation_transcription', False),
		variable_elongation_translation=args.get('variable_elongation_translation', False),
		disable_ribosome_capacity_fitting=not args.get('ribosome_fitting', True),
		disable_rnapoly_capacity_fitting=not args.get('rnapoly_fitting', True))

	task.run_task({})


class BaseSolver(abc.ABC):
	def __init__(self, method, args):
		self._method = method
		self._sim_dir = args.sim_path
		self._cpus = args.cpus
		self.perturbations = {}

		self.learning_rate = args.learning_rate
		self.parameter_step = args.parameter_step
		self.max_change = args.max_change  # TODO: implement generally

		self.iteration = args.starting_iteration
		self.variant = self.iteration * self.n_variants_per_iteration()

	@abc.abstractmethod
	def get_parameter_updates(self, original_values: Dict, objectives: List[float], paths: List[str]) -> Dict:
		"""
		Get the new parameter values based on the results of the current
		iteration.

		Args:
			original_values: the original values of parameters at the start of
				the iteration
			objectives: objective values for each variant run in this iteration
			paths: paths to the data objects (raw_data or sim_data) for each
				variant in this iteration

		Returns:
			mapping of each parameter to be modified to its update (the amount
			to change the parameter by)
		"""
		pass

	@abc.abstractmethod
	def get_parameter_perturbations(self, index: int) -> Tuple[Dict, Dict]:
		"""
		New parameter values for the given parameter set index within this
		iteration.

		Args:
			index: relative variant within the current iteration (0 is the first
				paramter set, 1 is the second and so on) so that parameters can
				be modified in different ways for each parameter set

		Returns:
			mapping of each parameter to be modified to its new value
		"""
		pass

	@abc.abstractmethod
	def n_variants_per_iteration(self) -> int:
		"""Number of variants (modified parameter sets) for each iteration"""
		pass

	def perturb_parameters(self, variants: List[int], raw_data_file: str, sim_data_file: str) -> List[str]:
		"""
		Perturb parameters for the upcoming iteration and run the parca
		if necessary (raw_data is modified).  Can run the parca in parallel
		for different sets of parameters.

		Args:
			variants: variant numbers for the current iteration to get the
				associated data files
			raw_data_file: path to the previous raw_data_file that can be used
				as a basis for updating to the new parameters
			sim_data_file: path to the previous sim_data_file that can be used
				as a basis for updating to the new parameters, not used
				if the parca needs to be rerun

		Returns:
			sim_data_files: paths to the newly created sim_data files for each
				variant
		"""

		sim_data_files = []
		sim_data_updates = []
		results = []

		pool = parallelization.pool(num_processes=self._cpus, nestable=True)
		cpus_per_parca = max(1, int(self._cpus / len(variants)))
		for variant in variants:
			index = variant - variants[0]
			raw_updates, sim_updates = self.get_parameter_perturbations(index)
			self.perturbations[(self.iteration, index)] = (raw_updates, sim_updates)
			new_raw_data_file, new_sim_data_file, metrics_file = self.data_paths(variant)
			if raw_updates:
				self.apply_updates(raw_data_file, raw_updates, new_raw_data_file)
				results.append(pool.apply_async(run_parca, (self._method.parca_args,
					new_raw_data_file, new_sim_data_file, metrics_file, cpus_per_parca)))
				sim_data_updates.append(sim_updates)
			else:
				self.apply_updates(sim_data_file, sim_updates, new_sim_data_file)
			sim_data_files.append(new_sim_data_file)

		pool.close()
		pool.join()
		pool = None

		for result, filename, updates in zip(results, sim_data_files, sim_data_updates):
			result.get()  # Get results in case an error was raised
			self.apply_updates(filename, updates, filename)

		return sim_data_files

	def update_parameters(self, variants: List[int], objectives: List[float]):
		"""
		Update parameters for the next iteration based on the objective results
		from the current iteration.

		Args:
			variants: variant numbers for the current iteration to get the
				associated data files
			objectives: objective values corresponding to each variant
		"""

		def update(data, objectives, paths):
			for param, update in self.get_parameter_updates(data, objectives, paths).items():
				original_value = data[param]
				if np.abs(units.strip_empty_units(update / original_value)) > self.max_change:
					update = np.sign(units.strip_empty_units(update / original_value)) * self.max_change * original_value
				data[param] += update

		raw_data_paths = []
		sim_data_paths = []
		for variant in variants:
			raw_data, sim_data, _ = self.data_paths(variant)
			raw_data_paths.append(raw_data)
			sim_data_paths.append(sim_data)

		update(self._method.raw_params, objectives, raw_data_paths)
		update(self._method.sim_params, objectives, sim_data_paths)

	def run_sims(self, sim_params: List[Dict]) -> List[str]:
		"""
		Run the set of simulations for a given iteration with the options
		specified in sim_params.  Will run a simulation for each entry in
		sim_params

		Args:
			sim_params: args for running a simulation

		Returns:
			sim_dirs: paths to the simulation directories that were run
		"""

		pool = parallelization.pool(self._cpus)
		results = [pool.apply_async(run_sim, (p,)) for p in sim_params]
		pool.close()
		pool.join()

		sim_dirs = [result.get() for result in results]

		return sim_dirs

	def run_iteration(self):
		"""
		Run one iteration of the algorithm.  Perturb parameters, run parca,
		run sims, get objectives and update to new parameters for the next
		iteration.
		"""

		variants = list(range(self.variant, self.variant+self.n_variants_per_iteration()))
		raw_data_file, sim_data_file = self.get_reference_parameters()
		if not self._method.initialized:
			self._method.initialize(raw_data_file, sim_data_file, self.iteration)
		sim_data_files = self.perturb_parameters(variants, raw_data_file, sim_data_file)
		sim_params = self._method.get_sim_params(self._sim_dir, variants)
		sim_out_dirs = self.run_sims(sim_params)
		objectives = self._method.get_objective(sim_out_dirs, sim_data_files)
		self.update_parameters(variants, objectives)

		self.iteration += 1
		self.variant = variants[-1] + 1

		return objectives

	def print_update(self, objectives):
		print(f'Objectives: {objectives}')
		self._method.print_update()

	def get_param(self, param: Parameter, path: str) -> Any:
		"""
		Get the parameter value from a specified data (raw_data or sim_data)
		file.

		Args:
			param: the param object that specifies how to access the desired data
			path: path to the data pickle to load

		Returns:
			parameter value stored in the data object
		"""

		with open(path, 'rb') as f:
			obj = pickle.load(f)

		return param.get_param(obj)

	def apply_updates(self, old_path: str, updates: Dict[Parameter, Any], new_path: str):
		"""
		Apply parameter updates to attributes in a pickle object and save the new object.

		Args:
			old_path: path to the old pickle object that will be modified
			updates: updates to apply to attributes in the old pickle object,
				nested attributes should be separated by '.'
			new_path: path to the new pickle object to store the modifications
		"""

		with open(old_path, 'rb') as f:
			obj = pickle.load(f)

		if updates:
			print('Updating values:')
			for param, val in updates.items():
				print(f'\t{param}: {val}')
				param.set_param(obj, val)

		with open(new_path, 'wb') as f:
			pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)

	def get_reference_parameters(self) -> Tuple[str, str]:
		"""
		Get paths to raw_data and sim_data files for the reference variant.

		Returns:
			raw_data: path to the reference raw_data pickle object for this variant
			sim_data: path to the reference sim_data pickle object for this variant
		"""

		if self.variant == 0:
			raw_data, sim_data, metrics = self.data_paths()
			if not os.path.exists(raw_data):
				InitRawDataTask(output=raw_data).run_task({})
			if not os.path.exists(sim_data):
				run_parca(self._method.parca_args, raw_data, sim_data, metrics, self._cpus)
		else:
			raw_data, sim_data, _ = self.data_paths(self.variant - 1)

		return raw_data, sim_data

	def data_paths(self, variant: Optional[int] = None) -> Tuple[str, str, str]:
		"""
		Get paths to raw_data and sim_data files for a given variant.

		Args:
			variant: variant number to get the path for, if variant is None,
				gets the default, unmodified files

		Returns:
			raw_data: path to the raw_data pickle object for this variant
			sim_data: path to the sim_data pickle object for this variant
			metrics: path to the output metrics pickle object for this variant
		"""

		if variant is None:
			kb_dir = fp.makedirs(self._sim_dir, constants.KB_DIR)
			sim_data_filename = constants.SERIALIZED_SIM_DATA_FILENAME
		else:
			kb_dir = fp.makedirs(self.get_variant_dir(variant), constants.VKB_DIR)
			sim_data_filename = constants.SERIALIZED_SIM_DATA_MODIFIED

		raw_data = os.path.join(kb_dir, constants.SERIALIZED_RAW_DATA)
		sim_data = os.path.join(kb_dir, sim_data_filename)
		metrics = os.path.join(kb_dir, constants.SERIALIZED_METRICS_DATA_FILENAME)

		return raw_data, sim_data, metrics

	def get_variant_dir(self, variant: int) -> str:
		return os.path.join(self._sim_dir, f'{self._method.variant_name}_{variant:06n}')
