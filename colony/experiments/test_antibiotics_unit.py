import numpy as np

from colony.experiments.antibiotics import (
    get_antibiotics_timeline,
    ANTIBIOTIC_KEY,
)

class TestGetTimeline:

	def assert_timelines_equal(self, timeline, expected_timeline):
		assert len(timeline) == len(expected_timeline)
		for (actual_time, actual), (expected_time, expected) in zip(
			timeline, expected_timeline
		):
			assert actual_time == expected_time
			if expected == {}:
				assert actual == expected
			else:
				assert (
					actual[('fields', ANTIBIOTIC_KEY)].tolist()
					== expected[('fields', ANTIBIOTIC_KEY)].tolist()
				)
			assert len(actual) == len(expected)

	def test_one_pulse_one_bin(self):
		timeline = get_antibiotics_timeline(
			(1, 1), (1, 1), [(1, 1, 1)], 5)
		expected_timeline = [
			(1, {('fields', ANTIBIOTIC_KEY): np.ones((1, 1))}),
			(2, {('fields', ANTIBIOTIC_KEY): np.zeros((1, 1))}),
			(5, {}),
		]
		assert timeline == expected_timeline

	def test_selects_correct_bin(self):
		timeline = get_antibiotics_timeline(
			(2, 2), (1, 1), [(1, 1, 1)], 5)
		expected_timeline = [
			(1, {('fields', ANTIBIOTIC_KEY): np.ones((2, 2))}),
			(2, {('fields', ANTIBIOTIC_KEY): np.zeros((2, 2))}),
			(5, {}),
		]
		self.assert_timelines_equal(timeline, expected_timeline)

	def test_two_pulses(self):
		timeline = get_antibiotics_timeline(
			(1, 1),
			(1, 1),
			[(1, 1, 1), (4, 2, 2)],
			10,
		)
		expected_timeline = [
			(1, {('fields', ANTIBIOTIC_KEY): np.ones((1, 1))}),
			(2, {('fields', ANTIBIOTIC_KEY): np.zeros((1, 1))}),
			(4, {('fields', ANTIBIOTIC_KEY): np.full((1, 1), 2)}),
			(6, {('fields', ANTIBIOTIC_KEY): np.zeros((1, 1))}),
			(10, {}),
		]
		self.assert_timelines_equal(timeline, expected_timeline)
