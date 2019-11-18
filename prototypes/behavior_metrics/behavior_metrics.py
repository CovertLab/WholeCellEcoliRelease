#!/usr/bin/env python2

"""Check that model behavior metrics are within expected bounds
"""

from __future__ import absolute_import, division, print_function
from os import path
from typing import Dict, Any, List

import numpy as np
import pandas as pd

from prototypes.behavior_metrics.tablereader_utils import read_subcolumn
from prototypes.behavior_metrics.dependency_graph import DependencyGraph
from wholecell.io.tablereader import TableReader
from wholecell.utils import filepath


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


#: Path from repository root to metrics configuration JSON file
METRICS_CONF_PATH = "prototypes/behavior_metrics/metrics.json"

#: Path to the simulation output directory from repository root
SIM_OUT_DIR = (
	"out/manual/wildtype_000000/000000/generation_000000/000000/simOut"
)

#: Map from mode names to the functions that handle the mode
MODE_FUNC_MAP = {
	"end_start_ratio": calc_end_start_ratio,
	"mean": np.mean,
	"stdev": np.std,
	"min": np.min,
	"max": np.max,
	"add_two_arrays": np.add,
	"calc_active_fraction": calc_active_fraction,
	"elementwise_min": np.minimum,
}


class BehaviorMetrics(object):
	"""Tests for model behavior metrics"""

	def __init__(self, metrics_conf_path, sim_out_dir):
		# type: (str, str) -> None
		"""Store metrics configuration path and simulation output dir

		Arguments:
			metrics_conf_path: Path to the metrics configuration JSON
				file
			sim_out_dir: Path to the simulation output directory
		"""
		self.metrics_conf_path = metrics_conf_path
		self.sim_out_dir = sim_out_dir

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
		metrics_conf = filepath.read_json_file(self.metrics_conf_path)
		results = []
		for metric, config in metrics_conf.items():
			data = self.load_data_from_config(config["data"])
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
					result = pd.DataFrame(
						[[metric, op, expected_min, expected_max, metric_val]],
						columns=[
							"metric", "mode", "expected_min",
							"expected_max", "value"
						],
					)
					results.append(result)
		results_df = pd.concat(results, ignore_index=True)
		results_df["pass"] = (
			(results_df["expected_min"] <= results_df["value"])
			& (results_df["value"] <= results_df["expected_max"])
		)
		return results_df

	@staticmethod
	def _calculate_operation(op_config, data):
		# type: (Dict[str, Any], Dict[str, Any]) -> Any
		op_func = MODE_FUNC_MAP[op_config["function"]]
		func_args = [
			data[arg_label] for arg_label in op_config["args"]
		]
		return op_func(*func_args)

	def load_data_from_config(self, data_conf_json):
		# type: (Dict[str, Any]) -> Dict[str, Any]
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
					"subcolumn": "CPLX0-3953[c]"
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
				},
				"B": {
					"table": "tRNA",
					"// If only column specified, there should be only"
					"// 1 subcolumn, which will be returned as a vector"
					"column": "counts"
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

		Arguments:
			data_conf_json: Parsed JSON dictionary that defines any number
				of data sources and how to load data from them.

		Returns:
			A dictionary where each key is the name of a data source and
			each value is the loaded data for the key.

		Raises:
			InvalidDependencyGraphError: If the dependency graph created by the
			operation attributes contains any cycles.
		"""
		loaded_data = {}
		for source_name, source_config in data_conf_json.items():
			if "subcolumn" in source_config:
				data = read_subcolumn(
					self.sim_out_dir, source_config["table"],
					source_config["column"], source_config["subcolumn"]
				)
				loaded_data[source_name] = data
				continue
			reader = TableReader(
				path.join(self.sim_out_dir, source_config["table"]))
			if "column" in source_config:
				data = reader.readColumn(source_config["column"])
				loaded_data[source_name] = data
			elif "attribute" in source_config:
				data = reader.readAttribute(source_config["attribute"])
				loaded_data[source_name] = data
			else:
				raise ValueError(
					"{} has neither 'column' nor 'attribute'".format(source_config))
		return loaded_data

	@staticmethod
	def order_operations(operation_configs):
		# type: (Dict[str, Dict[str, Any]]) -> List[str]
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
				deps = [
					arg for arg in config["args"]
					if arg in operation_configs.keys()
				]
				for dep in deps:
					graph.add_dep_relation(op_name, dep)
		return graph.get_topological_ordering()


def main():
	"""Main function that runs tests"""
	metrics = BehaviorMetrics(METRICS_CONF_PATH, SIM_OUT_DIR)
	results = metrics.calc_metrics()
	pd.options.display.width = None
	print(results)


if __name__ == "__main__":
	main()
