#!/usr/bin/env python3

"""Check that model behavior metrics are within expected bounds
"""

from __future__ import absolute_import, division, print_function
import json
from os import path

import numpy as np
import pandas as pd


from wholecell.io.tablereader import TableReader


"""Path from repository root to metrics configuration JSON file"""
METRICS_CONF_PATH = "prototypes/behavior_metrics/metrics.json"

"""Path to the simulation output directory from repository root"""
SIM_OUT_DIR = (
	"out/manual/wildtype_000000/000000/generation_000000/000000/simOut"
)


class BehaviorMetrics(object):
	"""Tests for model behavior metrics"""

	def __init__(self, metrics_conf_path, sim_out_dir):
		# type: (str, str)
		"""Load metrics configuration and define modes

		Modes are defined by a map from mode string to the function that
		handles the mode.

		Configuration file ile must be structured as follows (in JSON):
			{
				"metric_name": {
					"table": "A folder in SIM_OUT_DIR"
					"column": "A file in SIM_OUT_DIR/table/"
					"mode": "A key from TestBehaviorMetrics.MODE_FUNC_MAP"
					"range": [min, max]
				}
			}

		Arguments:
			metrics_conf_path: Path to the metrics configuration JSON
				file
			sim_out_dir: Path to the simulation output directory
		"""
		with open(metrics_conf_path, "r") as f:
			self.metrics_conf = json.load(f)
		self.MODE_FUNC_MAP = {
			"end_start_ratio": self._calc_end_start_ratio,
		}
		self.sim_out_dir = sim_out_dir

	def test_metrics(self):
		# type: () -> pd.DataFrame
		"""Test all behavior metrics

		For each metric defined in the configuration file, calls the
		appropriate mode function.
		"""
		results = []
		for metric, config in self.metrics_conf.items():
			reader = TableReader(path.join(self.sim_out_dir, config["table"]))
			data = reader.readColumn(config["column"])
			metric_val = self.MODE_FUNC_MAP[config["mode"]](data)
			expected_min, expected_max = config["range"]
			result = pd.DataFrame(
				[[metric, expected_min, expected_max, metric_val]],
				columns=["metric", "expected_min", "expected_max", "value"])
			results.append(result)
		results_df = pd.concat(results, ignore_index=True)
		results_df["pass"] = (
			(results_df["expected_min"] <= results_df["value"])
			& (results_df["value"] <= results_df["expected_max"])
		)
		return results_df

	def _calc_end_start_ratio(self, data):
		# type: (np.ndarray) -> float
		"""Check the ratio of a data's end to its start

		Arguments:
			data: A 1-dimensional array of timeseries data, in
				chronological order.

		Returns:
			The ratio end / start.
		"""
		return data[-1] / data[0]


def main():
	metrics = BehaviorMetrics(METRICS_CONF_PATH, SIM_OUT_DIR)
	results = metrics.test_metrics()
	print(results)


if __name__ == "__main__":
	main()
