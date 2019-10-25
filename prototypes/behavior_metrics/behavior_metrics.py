#!/usr/bin/env python3

"""Check that model behavior metrics are within expected bounds
"""

from __future__ import absolute_import, division, print_function
import json
from os import path
from typing import Dict, Any

import numpy as np
import pandas as pd


from wholecell.io.tablereader import TableReader


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


"""Path from repository root to metrics configuration JSON file"""
METRICS_CONF_PATH = "prototypes/behavior_metrics/metrics.json"

"""Path to the simulation output directory from repository root"""
SIM_OUT_DIR = (
	"out/manual/wildtype_000000/000000/generation_000000/000000/simOut"
)

"""Map from mode names to the functions that handle the mode"""
MODE_FUNC_MAP = {
	"end_start_ratio": calc_end_start_ratio,
	"mean": np.mean,
	"stdev": np.std,
	"min": np.min,
	"max": np.max,
}


class BehaviorMetrics(object):
	"""Tests for model behavior metrics"""

	def __init__(self, metrics_conf_path, sim_out_dir):
		# type: (str, str)
		"""Store metrics configuration path and simulation output dir

		Arguments:
			metrics_conf_path: Path to the metrics configuration JSON
				file
			sim_out_dir: Path to the simulation output directory
		"""
		self.metrics_conf_path = metrics_conf_path
		self.sim_out_dir = sim_out_dir

	def test_metrics(self):
		# type: () -> pd.DataFrame
		"""Test all behavior metrics

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
		with open(self.metrics_conf_path, "r") as f:
			metrics_conf = json.load(f)
		results = []
		for metric, config in metrics_conf.items():
			data = self.load_data_from_config(config["data"])
			for mode, mode_config in config["modes"].items():
				mode_func = MODE_FUNC_MAP[mode_config["function"]]
				func_args = [
					data[arg_label] for arg_label in mode_config["args"]
				]
				metric_val = mode_func(*func_args)
				expected_min, expected_max = mode_config["range"]
				result = pd.DataFrame(
					[[metric, mode, expected_min, expected_max, metric_val]],
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

	def load_data_from_config(self, data_conf_json):
		# type: (Dict[str, Any]) -> Dict[str, Any]
		"""Load data as specified in a configuration JSON.

		Arguments:
			data_conf_json: Parsed JSON dictionary that defines any number
				of data sources and how to load data from them.

		Returns:
			A dictionary where each key is the name of a data source and
			each value is the loaded data for the key.
		"""
		loaded_data = {}
		for source_name, source_config in data_conf_json.items():
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


def main():
	metrics = BehaviorMetrics(METRICS_CONF_PATH, SIM_OUT_DIR)
	results = metrics.test_metrics()
	pd.options.display.width = None
	print(results)


if __name__ == "__main__":
	main()
