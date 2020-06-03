'''
AnalysisPaths: object for easily accessing file paths to simulations based on
variants, seeds, and generation.
'''

from __future__ import absolute_import, division, print_function

from os import listdir
from os.path import isdir, join
from re import match, findall
from itertools import chain
from typing import List
import numpy as np

from wholecell.utils import constants

class AnalysisPaths(object):
	'''
	Object for easily accessing file paths to simulations based on
	variants, seeds, and generation. Within a specified variant you
	can then access all seeds and/or generations.

	Example:
		from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
		ap = AnalysisPaths(simOutDir)
		ap.get_cells(variant = [0,3], seed = [0], generation = [1,2])

	Above example should return all paths corresponding to variants
	0 and 3, seed 0, and generations 1 and 2. If a field is left blank
	it is assumed that all values for that field are desired. If all
	fields are left blank all cells will be returned.
	'''

	def __init__(self, out_dir, variant_plot = False, multi_gen_plot = False, cohort_plot = False):
		assert sum((variant_plot, multi_gen_plot, cohort_plot)) == 1, "Must specify exactly one plot type!"

		generation_dirs = []
		if variant_plot:
			# Final all variant files
			all_dirs = listdir(out_dir)
			variant_out_dirs = []
			# Consider only those directories which are variant directories
			for directory in all_dirs:
				# Accept directories which have a string, an underscore, and then a string
				# of digits exactly 6 units long
				if match(r'.*_\d{6}$', directory) is not None:
					variant_out_dirs.append(join(out_dir, directory))

			# Check to see if only wildtype variant exists
			if len(variant_out_dirs) == 0:
				for directory in all_dirs:
					if directory.startswith("wildtype_"):
						variant_out_dirs.append(join(out_dir, directory))

			if len(variant_out_dirs) == 0:
				raise Exception("Variant directory specified has no variant files in it!")

			# Get all seed directories in each variant directory
			seed_out_dirs = []
			for variant_dir in variant_out_dirs:
				all_dirs = listdir(variant_dir)
				for directory in all_dirs:
					if match(r'^\d{6}$', directory) is not None:
						seed_out_dirs.append(join(variant_dir, directory))

			# Get all generation files for each seed
			generation_dirs = []
			for seed_dir in seed_out_dirs:
				generation_dirs.extend(chain.from_iterable(self._get_generations(seed_dir)))

		elif cohort_plot:
			# Get all seed directories in each variant directory
			seed_out_dirs = []
			all_dirs = listdir(out_dir)
			for directory in all_dirs:
				if match(r'^\d{6}$', directory) is not None:
					seed_out_dirs.append(join(out_dir, directory))

			# Get all generation files for each seed
			generation_dirs = []
			for seed_dir in seed_out_dirs:
				generation_dirs.extend(chain.from_iterable(self._get_generations(seed_dir)))


		elif multi_gen_plot:
			generation_dirs = list(chain.from_iterable(self._get_generations(out_dir)))

		self._path_data = np.zeros(len(generation_dirs), dtype=[
			("path", "a500"),
			("variant", "i8"),
			("seed", "i8"),
			("generation", "i8"),
			("variantkb", "a500")
			])

		generations = []
		seeds = []
		variants = []
		variant_kb = []
		for filePath in generation_dirs:
			# Find generation
			matches = findall(r'generation_\d{6}', filePath)
			if len(matches) > 1:
				raise Exception("Expected only one match for generation!")
			generations.append(int(matches[0][-6:]))

			# Find seed
			seeds.append(int(filePath[filePath.rfind('generation_')-7:filePath.rfind('generation_')-1]))

			# Find variant
			variants.append(int(filePath[filePath.rfind('generation_')-14:filePath.rfind('generation_')-8]))

			# Find variant kb
			variant_kb.append(
				join(filePath[:filePath.rfind('generation_') - 8], "kb",
					constants.SERIALIZED_SIM_DATA_MODIFIED))

		self._path_data["path"] = generation_dirs
		self._path_data["variant"] = variants
		self._path_data["seed"] = seeds
		self._path_data["generation"] = generations
		self._path_data["variantkb"] = variant_kb

		self.n_generation = len(set(generations))
		self.n_variant = len(set(variants))
		self.n_seed = len(set(seeds))

	def get_cells(self, variant = None, seed = None, generation = None):
		"""Returns file paths for all the simulated cells."""
		# TODO: Rename this to get_cell_paths()?
		if variant is None:
			variantBool = np.ones(self._path_data.shape)
		else:
			variantBool = self._set_match("variant", variant)

		if seed is None:
			seedBool = np.ones(self._path_data.shape)
		else:
			seedBool = self._set_match("seed", seed)

		if generation is None:
			generationBool = np.ones(self._path_data.shape)
		else:
			generationBool = self._set_match("generation", generation)

		return self._path_data['path'][np.logical_and.reduce((variantBool, seedBool, generationBool))]

	def get_variant_kb(self, variant):
		kb_path = np.unique(self._path_data['variantkb'][np.where(self._path_data["variant"] == variant)])
		assert kb_path.size == 1
		return kb_path[0]

	def get_variants(self):
		return sorted(np.unique(self._path_data["variant"]))

	def _get_generations(self, directory):
		# type: (str) -> List[List[str]]
		generation_files = [
			join(directory, f) for f in listdir(directory)
			if isdir(join(directory, f)) and "generation" in f]  # type: List[str]
		generations = [[] for _ in generation_files]  # type: List[List[str]]
		for gen_file in generation_files:
			generations[int(gen_file[gen_file.rfind('_') + 1:])] = self._get_individuals(gen_file)
		return generations

	def _get_individuals(self, directory):
		# type: (str) -> List[str]
		individual_files = [
			join(directory, f) for f in listdir(directory)
			if isdir(join(directory, f))]  # type: List[str]
		individuals = [''] * len(individual_files)  # type: List[str]
		for ind_file in individual_files:
			individuals[int(ind_file[ind_file.rfind('/') + 1:])] = ind_file
		return individuals

	def _set_match(self, field, value):
		union = np.zeros(self._path_data[field].size)
		for x in value:
			union = np.logical_or(self._path_data[field] == x, union)
		return union
