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


def plot(experiment_id):
	with open (SECRETS_PATH, 'r') as f:
		secrets = json.load(f)
	emitter_config = get_atlas_database_emitter_config(
		**secrets['database'])
	uri = emitter_config['host']
	client = MongoClient(uri)
	query = {'experiment_id': experiment_id}

	config_collection = client[
		emitter_config['database']
	].configuration
	history_collection = client[
		emitter_config['database']
	].history

	environment_config = config_collection.find_one({
		'experiment_id': experiment_id,
		'type': 'environment_config',
	})

	unique_time_objs = history_collection.aggregate([
		{
			'$match': {
				'experiment_id': experiment_id
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

	# Plot snapshots
	agents = {
		timepoint['time']: timepoint['agents'] for timepoint in data}
	fields = {
		timepoint['time']: timepoint['fields'] for timepoint in data}
	snapshots_data = {
		'agents': agents,
		'fields': fields,
		'config': environment_config,
	}
	out_dir = os.path.join(OUT_DIR, experiment_id)
	os.makedirs(out_dir)
	plot_config = {
		'out_dir': out_dir,
		'filename': 'snapshot',
	}
	plot_snapshots(snapshots_data, plot_config)

	# Plot Emitted Values
	plot_settings = {
		'agents_key': 'agents',
	}


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument(
		"experiment_id",
		help="Experiment ID as recorded in the database"
	)
	args = parser.parse_args()
	plot(args.experiment_id)


if __name__ == '__main__':
	main()
