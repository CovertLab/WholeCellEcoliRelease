'''
AnalysisPaths: object for easily accessing file paths to simulations based on
variants, seeds, and generation.
'''

from __future__ import annotations

from os import listdir
from os.path import isdir, join
import re
from itertools import chain
from typing import cast, Iterable, List, Optional, Union

import numpy as np

from wholecell.utils import constants

class AnalysisPaths(object):
	'''
	Object for easily accessing file paths to simulations based on
	variants, seeds, and generation. Within a specified variant you
	can then access all seeds and/or generations.

	Example:
		from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
		ap = AnalysisPaths(simOutDir, variant_plot=True)
		ap.get_cells(variant = [0,3], seed = [0], generation = [1,2])

	Above example should return all paths corresponding to variants
	0 and 3, seed 0, and generations 1 and 2. If a field is left blank
	it is assumed that all values for that field are desired. If all
	fields are left blank all cells will be returned.

	* For a variant_plot, out_dir must be a top level simulation output dir.
	* For a cohort_plot, out_dir must be a variant output dir.
	* For a multi_gen_plot, out_dir must be a seed output dir.
	'''
	VARIANT_PATTERN = re.compile(r'.+_\d{6}')
	SEED_PATTERN = re.compile(r'\d{6}')

	def __init__(self, out_dir, *,
				 variant_plot: bool = False, multi_gen_plot: bool = False,
				 cohort_plot: bool = False) -> None:
		assert variant_plot + multi_gen_plot + cohort_plot == 1, (
			"Must specify exactly one plot type!")

		generation_dirs = []  # type: List[str]
		if variant_plot:
			# Find all variant directories in the given simulation output dir
			all_dirs = listdir(out_dir)
			variant_out_dirs = []
			for directory in all_dirs:
				if self.VARIANT_PATTERN.fullmatch(directory):
					variant_out_dirs.append(join(out_dir, directory))

			# Check to see if only wildtype variants exist that didn't match the pattern
			if len(variant_out_dirs) == 0:
				for directory in all_dirs:
					if directory.startswith("wildtype_"):
						variant_out_dirs.append(join(out_dir, directory))

			if len(variant_out_dirs) == 0:
				raise Exception("Variant out_dir doesn't contain variants!")

			# Get all seed directories in each variant directory
			seed_out_dirs = []
			for variant_dir in variant_out_dirs:
				all_dirs = listdir(variant_dir)
				for directory in all_dirs:
					if self.SEED_PATTERN.fullmatch(directory):
						seed_out_dirs.append(join(variant_dir, directory))

			# Get all generation files for each seed
			generation_dirs = []
			for seed_dir in seed_out_dirs:
				generation_dirs.extend(chain.from_iterable(self._get_generations(seed_dir)))

		elif cohort_plot:
			# Find all seed directories in the given variant directory
			seed_out_dirs = []
			all_dirs = listdir(out_dir)
			for directory in all_dirs:
				if self.SEED_PATTERN.fullmatch(directory):
					seed_out_dirs.append(join(out_dir, directory))

			# Get all generation files for each seed
			generation_dirs = []
			for seed_dir in seed_out_dirs:
				generation_dirs.extend(chain.from_iterable(self._get_generations(seed_dir)))


		elif multi_gen_plot:
			# Find all generation directories in the given seed directory
			generation_dirs = list(chain.from_iterable(self._get_generations(out_dir)))

		self._path_data = np.zeros(len(generation_dirs), dtype=[
			("path", "U500"),
			("variant", "i8"),
			("seed", "i8"),
			("generation", "i8"),
			("variantkb", "U500")
			])

		generations = []
		seeds = []
		variants = []
		variant_kb = []
		for filePath in generation_dirs:
			# Find generation
			matches = re.findall(r'generation_\d{6}', filePath)
			if len(matches) > 1:
				raise Exception("Expected only one match for generation!")
			generations.append(int(matches[0][-6:]))

			# TODO(jerry): Parse the path instead of using hardwired string offsets.
			gen_subdir_index = filePath.rfind('generation_')

			# Extract the seed index
			# Assumes: 6-digit seed index SSSSSS in 'SSSSSS/generation_...'
			seeds.append(int(filePath[gen_subdir_index - 7 : gen_subdir_index - 1]))

			# Extract the variant index
			# Assumes: 6-digit variant index VVVVVV in 'VARIANT-TYPE_VVVVVV/SSSSSS/generation_...'
			variants.append(int(filePath[gen_subdir_index - 14 : gen_subdir_index - 8]))

			# Find the variant kb pickle
			variant_kb.append(
				join(filePath[: gen_subdir_index - 8],
					 constants.VKB_DIR, constants.SERIALIZED_SIM_DATA_MODIFIED))

		self._path_data["path"] = generation_dirs
		self._path_data["variant"] = variants
		self._path_data["seed"] = seeds
		self._path_data["generation"] = generations
		self._path_data["variantkb"] = variant_kb

		self.n_generation = len(set(generations))
		self.n_variant = len(set(variants))
		self.n_seed = len(set(seeds))

	def get_cells(self, variant = None, seed = None, generation = None):
		# type: (Optional[Iterable[Union[int, str]]], Optional[Iterable[int]], Optional[Iterable[int]]) -> np.ndarray
		"""Returns file paths for all the simulated cells matching the given
		variant number, seed number, and generation number collections, where
		None => all.
		"""
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
		# type: (Union[int, str]) -> str
		kb_path = np.unique(self._path_data['variantkb'][np.where(self._path_data["variant"] == variant)])
		assert kb_path.size == 1
		return kb_path[0]

	def get_variants(self):
		# type: () -> List[Union[int, str]]
		"""Return all the variant indexes."""
		return sorted(np.unique(self._path_data["variant"]))

	def get_cell_variant(self, path: str) -> int:
		"""Return the variant index for the given get_cells() sim path."""
		return self._path_data['variant'][self._path_index(path)]

	def get_cell_seed(self, path: str) -> int:
		"""Return the seed for the given get_cells() sim path."""
		return self._path_data['seed'][self._path_index(path)]

	def get_cell_generation(self, path: str) -> int:
		"""Return the generation number for the given get_cells() sim path."""
		return self._path_data['generation'][self._path_index(path)]

	def get_cell_variant_kb(self, path: str) -> str:
		"""Return the variant kb path (simData_Modified.cPickle) for the given
		get_cells() sim path.
		"""
		return self._path_data['variantkb'][self._path_index(path)]

	def _path_index(self, path: str) -> int:
		"""Return the index into _path_data for the given get_cells() sim path."""
		indexes = np.where(self._path_data['path'] == path)[0]
		assert indexes.size == 1
		return indexes[0]

	def _get_generations(self, directory):
		# type: (str) -> List[List[str]]
		"""Get a sorted list of the directory's generation paths, each as a list
		of daughter cell paths.
		ASSUMES: directory contains "generation_000000" thru "generation_GGGGGG".
		"""
		generation_files = [
			join(directory, f) for f in listdir(directory)
			if isdir(join(directory, f)) and "generation" in f]  # type: List[str]
		generations = [[] for _ in generation_files]  # type: List[List[str]]
		for gen_file in generation_files:
			generations[int(gen_file[gen_file.rfind('_') + 1:])] = self._get_individuals(gen_file)
		return generations

	def _get_individuals(self, directory):
		# type: (str) -> List[str]
		"""Get a sorted list of the directory's daughter cell paths, each of
		them a place for simOut/ and plotOut/ subdirs.
		ASSUMES: directory is a generation directory like "generation_000001",
		each containing numbered daughter subdirs "000000" thru "DDDDDD".
		"""
		individual_files = [
			join(directory, f) for f in listdir(directory)
			if isdir(join(directory, f))]  # type: List[str]
		individuals = [''] * len(individual_files)  # type: List[str]
		for ind_file in individual_files:
			individuals[int(ind_file[ind_file.rfind('/') + 1:])] = ind_file
		return individuals

	def _set_match(self, field, value):
		# type: (str, Iterable[Union[int, str]]) -> np.ndarray
		union = np.zeros(self._path_data[field].size)
		for x in value:
			union = cast(np.ndarray, np.logical_or(self._path_data[field] == x, union))
		return union
