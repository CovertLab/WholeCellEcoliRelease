'''Validate the data from an experiment.

Currently, we only check that rifampicin is conserved. For usage
information, run:
	python validator.py -h
'''

import argparse

import numpy as np

from vivarium.analysis.analyze import Analyzer
from vivarium.core.experiment import get_in
from vivarium.library.lattice_utils import get_bin_volume
from vivarium.library.units import units


DIMENSIONS_PATH = ('dimensions',)
FIELDS_PATH = ('fields',)
AGENTS_PATH = ('agents',)
AGENT_VOLUME_PATH = ('boundary', 'volume')
RIFAMPICIN_PATH = ('boundary', 'cytoplasm', 'rifampicin')


def assert_rifampicin_conserved(data):
	'''Check that the total moles of rifampicin is constant over time.

	At every timepoint, sum the moles of rifampicin in the environment
	and in all the cells. Assert that the minimum and maximum value are
	close (using np.allclose).
	'''
	moles = [
		get_timepoint_total_moles(RIFAMPICIN_PATH, data[time])
		for time in sorted(data.keys())
	]
	# Check that we successfully downloaded data
	assert len(moles) > 1, 'Fewer than 2 timepoints retrieved.'
	min_moles = min(moles)
	max_moles = max(moles)
	assert np.allclose(min_moles, max_moles), (
		'Moles of rifampicin not constant: {}'.format(moles))


def get_timepoint_total_moles(path_from_agent, timepoint_data):
	'''Get the total moles of a molecule present at a timepoint

	Sums the moles present in all the agents and in the environment.

	Arguments:
		path_from_agent (tuple): Path from the agent root to the
			variable holding the molecule's concentration. We assume
			that the same key identifies the molecule in the agent (i.e.
			the last element of this path) and in the environment
			fields.
		timepoint_data (dict): Data for a timepoint.

	Returns:
		Quantity: Number of molecules present, in units of moles.
	'''
	key = path_from_agent[-1]
	# Get environment dimensions
	dimensions = get_in(timepoint_data, DIMENSIONS_PATH)
	n_bins = dimensions['n_bins']
	bounds = dimensions['bounds']
	depth = dimensions['depth']
	bin_volume = get_bin_volume(n_bins, bounds, depth) * units.L

	# Sum counts across fields
	field = get_in(timepoint_data, FIELDS_PATH + (key,))
	conc_sum = np.sum(field) * units.mmol / units.L
	mol_environment = conc_sum * bin_volume

	# Sum counts across agents
	agents_data = get_in(timepoint_data, AGENTS_PATH)
	mol_agents = 0
	for _, agent_data in agents_data.items():
		concentration = get_in(agent_data, path_from_agent) * (
			units.mmol / units.L)
		agent_volume = get_in(agent_data, AGENT_VOLUME_PATH) * (
			units.fL)
		mol_agents += concentration * agent_volume

	# Return total moles
	return (mol_environment + mol_agents).to(units.mol)


def main():
	'''Validate data.'''
	parser = argparse.ArgumentParser()
	parser.add_argument(
		'experiment_id',
		type=str,
		help='ID of experiment to validate data from.'
	)
	Analyzer.add_connection_args(parser)
	args = parser.parse_args()

	data, _ = Analyzer.get_data(args, args.experiment_id)
	assert_rifampicin_conserved(data)


if __name__ == '__main__':
	main()
