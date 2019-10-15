#!/usr/bin/env python3

"""Check that model behavior metrics are within expected bounds
"""

from __future__ import absolute_import, division, print_function
import json
from os import path
import unittest

from wholecell.io.tablereader import TableReader


"""Path from repository root to metrics configuration JSON file

File must be structured as follows:
	{
		"metric_name": {
			"table": "A folder in SIM_OUT_DIR"
			"column": "A file in SIM_OUT_DIR/table/"
			"mode": "A key from TestBehaviorMetrics.MODE_FUNC_MAP"
			"range": [min, max]
		}
	}
"""
METRICS_CONF_PATH = "prototypes/behavior_metrics/metrics.json"

"""Path to the simulation output directory from repository root"""
SIM_OUT_DIR = (
	"out/manual/wildtype_000000/000000/generation_000000/000000/simOut"
)


class TestBehaviorMetrics(unittest.TestCase):
	"""Tests for model behavior metrics"""

	longMessage = True

	@classmethod
	def setUpClass(cls):
		"""Load metrics configuration and define modes

		Modes are defined by a map from mode string to the function that
		handles the mode.
		"""
		with open(METRICS_CONF_PATH, "r") as f:
			cls.metrics_conf = json.load(f)
		cls.MODE_FUNC_MAP = {
			"end_start_ratio": cls._check_end_start_ratio,
		}

	def test_metrics(self):
		"""Test all behavior metrics

		For each metric defined in the configuration file, calls the
		appropriate mode function.
		"""
		for metric, config in self.metrics_conf.items():
			reader = TableReader(path.join(SIM_OUT_DIR, config["table"]))
			data = reader.readColumn(config["column"])
			msg = "Metric: {}".format(metric)
			self.MODE_FUNC_MAP[config["mode"]](
				self, data, config["range"], msg)

	def _check_end_start_ratio(self, data, expected_range, msg):
		# (TestBehaviorMetrics, np.ndarray, Tuple[float, float]) -> None
		"""Check the ratio of a data's end to its start

		The ratio is computed as: end / start

		Arguments:
			data: A 1-dimensional array of timeseries data, in
				chronological order.
			expected_range: A tuple (min, max) that defines the minimum
				and maximum allowable values for the ratio.
			msg: String to display along with any test failures. Should
				at minimum include the metric being tested.
		"""
		end_start_ratio = data[-1] / data[0]
		self.assertLessEqual(expected_range[0], end_start_ratio, msg)
		self.assertLessEqual(end_start_ratio, expected_range[1], msg)


if __name__ == "__main__":
	unittest.main()
