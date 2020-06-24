from __future__ import absolute_import, division, print_function

import argparse
import json
import os

from pymongo import MongoClient
from vivarium.plots.multibody_physics import plot_snapshots
from vivarium.core.composition import plot_agents_multigen

from lattice_experiment import get_atlas_database_emitter_config
from lattice_experiment import SECRETS_PATH


OUT_DIR = 'out'


def plot(args):
	# Connect to mongoDB Atlas Database
	with open(SECRETS_PATH, 'r') as f:
		secrets = json.load(f)
	emitter_config = get_atlas_database_emitter_config(
		**secrets['database'])
	uri = emitter_config['host']
	client = MongoClient(uri)

	# Retrieve experiment data
	config_collection = client[
		emitter_config['database']
	].configuration
	history_collection = client[
		emitter_config['database']
	].history

	environment_config = config_collection.find_one({
		'experiment_id': args.experiment_id,
		'type': 'environment_config',
	})

	unique_time_objs = history_collection.aggregate([
		{
			'$match': {
				'experiment_id': args.experiment_id
			}
		}, {
			'$group': {
				'_id': {
					'time': '$time'
				},
				'id': {
					'$first': '$_id'
				}
			}
		}, {
			'$sort': {
				'_id.time': 1
			}
		},
	])
	unique_time_ids = [
		obj['id'] for obj in unique_time_objs
	]
	data_cursor = history_collection.find({
		'_id': {
			'$in': unique_time_ids
		}
	}).sort('time')
	data = list(data_cursor)

	# Create plots
	out_dir = os.path.join(OUT_DIR, args.experiment_id)
	if os.path.exists(out_dir):
		if not args.force:
			raise IOError('Directory {} already exists'.format(out_dir))
	else:
		os.makedirs(out_dir)

	if args.snapshots:
		agents = {
			timepoint['time']: timepoint['agents'] for timepoint in data
			if timepoint['time'] != 0
		}
		fields = {
			timepoint['time']: timepoint['fields'] for timepoint in data
			if timepoint['time'] != 0
		}
		snapshots_data = {
			'agents': agents,
			'fields': fields,
			'config': environment_config,
		}
		plot_config = {
			'out_dir': out_dir,
			'filename': 'snapshot',
		}
		plot_snapshots(snapshots_data, plot_config)

	if args.timeseries:
		plot_settings = {
			'agents_key': 'agents',
			'title_size': 10,
		}
		dict_data = {
			timepoint['time']: timepoint for timepoint in data
			if timepoint['time'] != 0
		}
		plot_agents_multigen(dict_data, plot_settings, out_dir)


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument(
		'experiment_id',
		help='Experiment ID as recorded in the database',
	)
	parser.add_argument(
		'--snapshots', '-s',
		action='store_true',
		default=False,
		help='Plot snapshots',
	)
	parser.add_argument(
		'--timeseries', '-t',
		action='store_true',
		default=False,
		help='Generate line plot for each variable over time',
	)
	parser.add_argument(
		'--force', '-f',
		action='store_true',
		default=False,
		help=(
			'Write plots even if output directory already exists. This '
			+ 'could overwrite your existing plots'
		),
	)
	args = parser.parse_args()
	plot(args)


if __name__ == '__main__':
	main()
