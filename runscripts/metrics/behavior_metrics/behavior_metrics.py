#!/usr/bin/env python2

"""Compute metrics for model behavior and check expected bounds
"""

from __future__ import absolute_import, division, print_function

from collections import namedtuple
import importlib
from os import path
import re
from typing import Any, Dict, Iterable, List, Optional, Sequence, Text, Tuple, Union

import numpy as np
import pandas as pd
from unum import Unum
from six.moves import cPickle, range

from wholecell.io.tablereader import TableReader
from wholecell.utils import filepath, units, toya
from wholecell.utils.dependency_graph import DependencyGraph
from wholecell.utils.protein_counts import (
	get_simulated_validation_counts,
)


# A type alias for Python 2 str or unicode; Python 3 str (not bytes).
# After porting to Python 3, we can use plain `str`.
# Note: Unions don't work with `isinstance()`.
String = Union[str, Text]


def calc_end_start_ratio(data):
	# type: (np.ndarray) -> float
	"""Check the ratio of a data's end to its start

	Arguments:
		data: A 1-dimensional array of timeseries data, in
			chronological order.

	Returns:
		The ratio end / start.
	"""
	return data[-1] / data[0]


def calc_active_fraction(active_counts, inactive_counts):
	# type: (np.ndarray, np.ndarray) -> np.ndarray
	"""Compute active fraction.

	Arguments:
		active_counts: Counts of active units over time.
		inactive_counts: Counts of inactive units at the same times as
			active_counts.

	Returns:
		At each time point, the active fraction.
	"""
	return active_counts / (active_counts + inactive_counts)


def find_limiting_metabolites(counts, names, window):
	# type: (np.ndarray, Sequence[str], int) -> Iterable[str]
	"""Find all metabolites that are limiting for some period of time.

	A metabolite is considered limiting over a window of time if, within
	that window, its count remains constant.

	Arguments:
		counts: Counts of the metabolites over time. Should be a matrix
			with a column for each metabolite and a row for each time
			point.
		names: Names of the metabolites, in the same order as the
			columns for the metabolites appear in counts.
		window: Length of the window to use, in units of rows of counts.

	Returns:
		An iterable collection with the unsorted names of all
		metabolites ever found to be limiting in counts.
	"""
	limiting = set()
	diff = np.diff(counts, axis=0)
	for i in range(diff.shape[0] - window):
		production_in_window = np.any(diff[i:i + window] > 0, axis=0)
		i_unproduced_metabolites = np.where(
			production_in_window == False)[0].astype(int)
		if len(i_unproduced_metabolites):
			curr_limiting = names[i_unproduced_metabolites]
			limiting.update(curr_limiting)
	return limiting


def normalize_to_column(to_normalize, column_index):
	# type: (np.ndarray, int) -> np.ndarray
	"""Normalize a matrix by dividing all columns by some column.

	Arguments:
		to_normalize: The matrix to normalize.
		column_index: The index of the column in the matrix to divide
			by.

	Returns:
		The normalized matrix.
	"""
	return to_normalize / to_normalize[column_index, :]


def find_indices_bulk(to_search, targets):
	# type: (List[Any], Iterable[Any]) -> np.ndarray
	"""Generate a list of the indices at which elements appear in list.

	Only the first occurrence of each element is included.

	Arguments:
		to_search: List to lookup elements in
		targets: The elements to lookup

	Returns:
		The indices.
	"""
	return np.array([to_search.index(elem) for elem in targets], np.int)


def np_pick(array, pick_spec):
	# type: (np.ndarray, Iterable[Union[List[int], int]]) -> Any
	"""Perform numpy indexing based on a specification

	Arguments:
		array: The array to slice.
		pick_spec: The selection specification, which should consist
			of a list of axis specifications. Each axis specification can
			either be an integer, to represent an index, or as a list of
			indices.

	Returns:
		The given array indexed by the specified indices.
	"""
	parsed_spec = []  # type: List[Union[Tuple[int, ...], int]]
	for axis_spec in pick_spec:
		if isinstance(axis_spec, list):
			parsed_spec.append(tuple(axis_spec))
		else:
			parsed_spec.append(axis_spec)
	return array[tuple(parsed_spec)]


#: Path from repository root to metrics configuration JSON file
METRICS_CONF_PATH = "runscripts/metrics/behavior_metrics/metrics.json"

#: Path to the simulation output directory from repository root
SIM_OUT_DIR = (
	"out/manual/wildtype_000000/000000/generation_000000/000000/simOut"
)

#: Path to pickle that stores validation data for correlation checks
VALIDATION_PICKLE_PATH = "out/manual/kb/validationData.cPickle"

#: Path to pickle that stores pickle of data for metrics computation
METRICS_PICKLE_PATH = "out/manual/kb/metricsData.cPickle"

#: Map from mode names to the functions that handle the mode
MODE_FUNC_MAP = {
	# Simple Scalar-Valued Functions
	"mean": np.mean,
	"stdev": np.std,
	"min": np.min,
	"max": np.max,
	"len": len,
	"sum": np.sum,
	"absolute": np.absolute,

	# Simple Binary Operators
	"*": lambda x, y: x * y,
	"/": lambda x, y: x / y,
	"-": lambda x, y: x - y,

	# Array Operations
	"add_two_arrays": np.add,
	"elementwise_min": np.minimum,
	"elementwise_divide": np.divide,
	"pairwise_diffs": np.diff,
	"pairwise_diffs_axis": lambda a, axis: np.diff(a, axis=axis),
	"slice": lambda arr, start, end: arr[start:end],
	"np_pick": np_pick,
	"last_elem": lambda x: x[-1],
	"ravel": np.ravel,

	# Metrics
	"end_start_ratio": calc_end_start_ratio,
	"calc_active_fraction": calc_active_fraction,
	"pearson_correlation": lambda x, y: np.corrcoef(x, y)[0, 1],
	"pearson_correlation_log10": lambda x, y: np.corrcoef(
		np.log10(x + 1), np.log10(y + 1))[0, 1],
	"find_limiting_metabolites": find_limiting_metabolites,

	# Processing
	"adjust_toya_data": toya.adjust_toya_data,
	"process_simulated_fluxes": (
		lambda filter_ids, rxn_ids, fluxes, id_map:
		toya.process_simulated_fluxes(
			filter_ids, rxn_ids, fluxes, id_map
		)[0]
	),
	"process_toya_data": toya.process_toya_data,
	"fluxome_common_ids": toya.get_common_ids,
	"fluxome_root_to_id_indices_map":
		toya.get_root_to_id_indices_map,
	"normalize_to_column": normalize_to_column,
	"get_simulated_validation_counts": get_simulated_validation_counts,
	"filter_no_nan": lambda a: a[~np.isnan(a)],
	"find_indices_bulk": find_indices_bulk,
	"strip_units": lambda to_strip, unit: to_strip.asNumber(unit),
}


class BehaviorMetrics(object):
	"""Tests for model behavior metrics"""

	def __init__(
		self, metrics_conf_path, sim_out_dir, validation_path=None,
		metrics_pickle_path=None
	):
		# type: (str, str, str, str) -> None
		"""Store provided paths.

		Arguments:
			metrics_conf_path: Path to the metrics configuration JSON
				file.
			sim_out_dir: Path to the simulation output directory.
			validation_path: Path to validation data cPickle. May be
				excluded or set to None if no validation data will be
				used.
			metrics_pickle_path: Path to pickle of data for metrics
				computation.
		"""
		self.metrics_conf_path = metrics_conf_path
		self.sim_out_dir = sim_out_dir
		self.validation_path = validation_path
		self.metrics_pickle_path = metrics_pickle_path

	def calc_metrics(self):
		# type: () -> pd.DataFrame
		"""Compute all behavior metrics and determine whether in bounds

		For each metric defined in the configuration file, calls the
		appropriate mode functions and checks that the results are
		within expected ranges.

		Returns:
			A table of the results with columns:
				* metric
				* mode
				* expected_min
				* expected_max
				* value
				* pass
			where metric and mode specify the test run. value is the
			result of the mode computation, and value is expected to be
			within [expected_min, expected_max]. Pass is a boolean
			representation of whether this is true.
		"""
		Result = namedtuple(
			"Result",
			[
				"metric", "mode", "expected_min", "expected_max",
				"value", "passes"
			],
		)
		metrics_conf = filepath.read_json_file(self.metrics_conf_path)
		pickles = {}
		if self.validation_path:
			with open(self.validation_path, "rb") as f:
				pickles["validation_data"] = cPickle.load(f)
		if self.metrics_pickle_path:
			with open(self.metrics_pickle_path, "rb") as f:
				pickles["metrics_data"] = cPickle.load(f)
		results = []
		for metric, config in metrics_conf.items():
			data = self.load_data_from_config(config["data"], pickles)
			ordered_ops = BehaviorMetrics.order_operations(
				config["operations"]
			)
			for op in ordered_ops:
				op_config = config["operations"][op]
				metric_val = BehaviorMetrics._calculate_operation(
					op_config, data)
				data[op] = metric_val
				if "range" in op_config:
					expected_min, expected_max = op_config["range"]
					in_bounds = (
						expected_min <= metric_val <= expected_max)
					result = Result(
						metric, op, expected_min, expected_max,
						metric_val, in_bounds
					)
					results.append(result)
				elif "expected_set" in op_config:
					expected = set(op_config["expected_set"])
					matches = expected == metric_val
					result = Result(
						metric, op, None, None,
						metric_val, matches
					)
					results.append(result)

		results_df = pd.DataFrame(results)
		return results_df

	@staticmethod
	def _resolve_func_arg(arg, data):
		if isinstance(arg, list):
			return [
				BehaviorMetrics._resolve_func_arg(elem, data) for elem in arg]
		return data[arg]

	@staticmethod
	def _calculate_operation(op_config, data):
		# type: (Dict[str, Any], Dict[String, Any]) -> Any
		op_func = MODE_FUNC_MAP[op_config["function"]]
		func_args = [
			BehaviorMetrics._resolve_func_arg(arg, data) for arg in op_config["args"]
		]
		return op_func(*func_args)

	def load_data_from_config(self, data_conf_json, pickles=None):
		# type: (Dict[String, Any], Optional[dict]) -> Dict[String, Any]
		"""Load data as specified in a configuration JSON.

		The configuration JSON should be structured as follows:

		"ribosomeActiveFraction": {
			"data": {
				"// active_counts is an arbitrary label used to specify"
				"// arguments to operations"
				"active_counts": {
					"table": "UniqueMoleculeCounts",
					"column": "uniqueMoleculeCounts",
					"subcolumn": "activeRibosome"
				},
				"s30_counts": {
					"table": "BulkMolecules",
					"column": "counts",
					"subcolumn": "CPLX0-3953[c]",
				},
				"s50_counts": {
					"table": "BulkMolecules",
					"column": "counts",
					"subcolumn": "CPLX0-3962[c]"
				},
				"A": {
					"table": "BulkMolecules",
					"// We can also read in attributes"
					"attribute": "subcolumns"
					"// Units can be specified as strings of attributes"
					"// of utils.units"
					"units": "g*m/L"
				},
				"B": {
					"table": "tRNA",
					"// If only column specified, there should be only"
					"// 1 subcolumn, which will be returned as a vector"
					"column": "counts"
					"// Units can also be specified by a name to import"
					"units": {
						"import": "model.ecoli.constants.DUMMY_UNIT"
					}
				}
			},
			"operations": {
				"// Operations can take any other operation as input (no cycles)"
				"inactive_counts": {
					"function": "sum",
					"args": ["s30_counts", "s50_counts"]
				},
				"active_fraction": {
					"function": "calc_active_fraction",
					"args": ["active_counts", "inactive_counts"]
				},
				"mean": {
					"function": "mean",
					"args": ["active_fraction"],
					"// Any operation with a range will have its output validated"
					"range": [0, 0]
				}
			}
		}

		Note that units are left-multiplied, so when applied to numpy
		vectors, they will apply to the whole vector, not each element.

		Arguments:
			data_conf_json: Parsed JSON dictionary that defines any number
				of data sources and how to load data from them.
			pickles: Dictionary of loaded pickles. Each key-value pair
				maps the name by which the pickle is referenced in the
				config to the loaded pickle.

		Returns:
			A dictionary where each key is the name of a data source and
			each value is the loaded data for the key.

		Raises:
			InvalidDependencyGraphError: If the dependency graph created by the
			operation attributes contains any cycles.
		"""
		if pickles is None:
			pickles = {}
		loaded_data = {}
		for source_name, source_config in data_conf_json.items():
			if "constant" in source_config:
				data = source_config["constant"]
			elif "subcolumn" in source_config:
				reader = TableReader(
					path.join(self.sim_out_dir, source_config["table"]))
				data = reader.readSubcolumn(
					source_config["column"], source_config["subcolumn"])
			elif "import" in source_config:
				data = BehaviorMetrics._load_from_import_string(
					source_config["import"])
			elif "cPickle" in source_config:
				pickle = pickles[source_config["cPickle"]]
				data = pickle
				if "dotted_name" in source_config:
					data = BehaviorMetrics._resolve_dotted_name(
						data, source_config["dotted_name"])
				if "dict_key" in source_config:
					data = data[source_config["dict_key"]]
			elif "table" in source_config:
				reader = TableReader(
					path.join(self.sim_out_dir, source_config["table"]))
				if "column" in source_config:
					data = reader.readColumn(source_config["column"])
				elif "attribute" in source_config:
					data = reader.readAttribute(source_config["attribute"])
				else:
					raise ValueError(
						"{} has neither 'column' nor 'attribute'".format(source_config))
			else:
				raise ValueError(
					"{} has none of 'constant', 'subcolumn', and 'table'".format(
						source_config))
			if "units" in source_config:
				parsed_units = BehaviorMetrics.parse_units(
					source_config["units"])
				data = parsed_units * data
			loaded_data[source_name] = data

		return loaded_data

	@staticmethod
	def _resolve_dotted_name(obj, name):
		# type: (object, String) -> Any
		names = name.split(".")
		for name_part in names:
			obj = getattr(obj, name_part)
		return obj

	@staticmethod
	def _flatten(lst):
		# type: (List[Any]) -> List[Any]
		"""Flatten a list"""
		flat = []
		for elem in lst:
			if isinstance(elem, list):
				flat.extend(BehaviorMetrics._flatten(elem))
			else:
				flat.append(elem)
		return flat

	@staticmethod
	def parse_units(unit_def):
		# type: (Union[String, Dict[String, Any]]) -> Unum
		"""Get an Unum object that can store the specified units

		Arguments:
			unit_def: A definition of a unit, either as a string or as a
				dict. Strings are processed by get_units_str, and dicts
				are processed by get_units_dict.

		Returns:
			An Unum object storing the specified units.
		"""
		if isinstance(unit_def, (str, Text)):
			return BehaviorMetrics.parse_units_str(unit_def)
		return BehaviorMetrics.parse_units_dict(unit_def)

	@staticmethod
	def parse_units_str(unit_str):
		# type: (String) -> Unum
		"""Get an Unum object from a unit string.

		Arguments:
			unit_str: A string specifying the unit. It is
				interpreted as one or more names of attributes of
				wholecell.utils.units, which may be combined only by
				multiplication (*) or division (/). No parentheses are
				allowed. The numeral 1 can be used to represent a
				unitless value, e.g. in 1/g.

		Returns:
			An Unum object storing the specified units.
		"""
		operator_regex = "[*/]"
		operators = re.findall(operator_regex, unit_str)
		atoms = re.split(operator_regex, unit_str)
		operations = zip(atoms[1:], operators)
		total_units = BehaviorMetrics._eval_atomic_unit_str(atoms[0])
		for atom, operator in operations:
			unit = BehaviorMetrics._eval_atomic_unit_str(atom)
			if operator == "*":
				total_units *= unit
			else:
				total_units /= unit
		return total_units

	@staticmethod
	def parse_units_dict(unit_dict):
		# type: (Dict[String, Any]) -> Unum
		"""Get an Unum object from a JSON object.

		Arguments:
			unit_dict: Dictionary representing a JSON object that
				defines how to retrieve the unit. For now, the only
				supported definition specifies a name to import using
				the `import` key.

		Returns:
			An Unum object storing the specified units.
		"""
		if "import" in unit_dict:
			return BehaviorMetrics._load_from_import_string(
				unit_dict["import"])
		else:
			raise ValueError(
				"The units specification {} was not recognized".format(
					unit_dict))

	@staticmethod
	def _load_from_import_string(import_str):
		"""Load a value from a name to import.

		The parent of the specified object must be a module importable
		by importlib.import_module.

		Arguments:
			import_str: Name to import.

		Returns:
			The imported value.
		"""
		split_name = import_str.split(".")
		parent_full_name = ".".join(split_name[:-1])
		parent = importlib.import_module(parent_full_name)
		return getattr(parent, split_name[-1])

	@staticmethod
	def _eval_atomic_unit_str(unit_str):
		# type: (String) -> Unum
		return 1 if unit_str == "1" else getattr(units, unit_str)

	@staticmethod
	def order_operations(operation_configs):
		# type: (Dict[String, Dict[String, Any]]) -> Sequence[String]
		"""Sorts operation configs for evaluation.

		Operations can take the results of other operations as input, so
		operations must be evaluated in an order such that when any
		operation is evaluated, its arguments have already been
		evaluated.

		Arguments:
			operation_configs: Dictionary with operation names as keys
				operation configurations as values.

		Returns:
			List of operation names in order of evaluation.

		Raises:
			InvalidDependencyGraphError: If the dependency graph created by the
			operation attributes contains any cycles.
		"""
		graph = DependencyGraph()
		graph.add_nodes(operation_configs.keys())
		for op_name, config in operation_configs.items():
			if "args" in config:
				args = BehaviorMetrics._flatten(config["args"])
				deps = [
					arg for arg in args
					if arg in operation_configs.keys()
				]
				for dep in deps:
					graph.add_dep_relation(op_name, dep)
		return graph.get_topological_ordering()


def main():
	"""Main function that runs tests"""
	metrics = BehaviorMetrics(
		METRICS_CONF_PATH, SIM_OUT_DIR, VALIDATION_PICKLE_PATH,
		METRICS_PICKLE_PATH)
	results = metrics.calc_metrics()
	pd.options.display.width = None
	print(results)


if __name__ == "__main__":
	main()
