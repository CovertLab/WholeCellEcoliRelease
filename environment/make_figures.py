'''Generate figures for wcEcoli colony simulation

For usage information, run:
	python make_figures.py -h
'''
import argparse
import os
import sys

from vivarium.analysis.analyze import Analyzer
from vivarium.plots.multibody_physics import plot_tags, plot_snapshots

import wholecell.utils.filepath as fp


TAG_PATH_NAME_MAP = {
	('boundary', 'bulk_molecules_report', 'XAPB-MONOMER[i]'): 'XapB',
	(
		'boundary', 'bulk_molecules_report', 'TRANS-CPLX-201[s]'
	): 'AcrAB-TolC',
}
OUT_DIR = os.path.join('environment', 'figs')
FILE_EXTENSION = 'pdf'
FIG_2_EXPERIMENT_ID = '20200806.170126'
FIG_3_EXPERIMENT_ID = '20200806.170126'
METADATA_FILE = 'metadata.json'


def get_metadata():
	metadata = {
		'git_hash': fp.run_cmdline('git rev-parse HEAD'),
		'git_branch': fp.run_cmdline('git symbolic-ref --short HEAD'),
		'time': fp.timestamp(),
		'python': sys.version.splitlines()[0],
		'fig2_id': FIG_2_EXPERIMENT_ID,
		'fig3_id': FIG_3_EXPERIMENT_ID,
	}
	return metadata


def make_fig2(data, environment_config):
	'''Figure shows heterogeneous expression within wcEcoli agents.'''
	tags_data = Analyzer.format_data_for_tags(data, environment_config)
	plot_config = {
		'out_dir': OUT_DIR,
		'tagged_molecules': TAG_PATH_NAME_MAP.keys(),
		'filename': 'fig2.{}'.format(FILE_EXTENSION),
		'tag_path_name_map': TAG_PATH_NAME_MAP,
		'tag_label_size': 36,
	}
	plot_tags(tags_data, plot_config)


def make_fig3(data, environment_config):
	'''Figure shows heterogeneous GLC uptake by across environment.'''
	snapshots_data = Analyzer.format_data_for_snapshots(
		data, environment_config)
	plot_config = {
		'out_dir': OUT_DIR,
		'filename': 'fig3.{}'.format(FILE_EXTENSION),
		'include_fields': ['rifampicin'],
		'field_label_size': 36,
	}
	plot_snapshots(snapshots_data, plot_config)


def main():
	'''Generate all figures.'''
	if not os.path.exists(OUT_DIR):
		os.makedirs(OUT_DIR)
	fp.write_json_file(os.path.join(
		OUT_DIR, METADATA_FILE), get_metadata())
	parser = argparse.ArgumentParser()
	Analyzer.add_connection_args(parser)
	args = parser.parse_args()

	data, environment_config = Analyzer.get_data(
		args, FIG_2_EXPERIMENT_ID)
	make_fig2(data, environment_config)

	if FIG_2_EXPERIMENT_ID != FIG_3_EXPERIMENT_ID:
		data, environment_config = Analyzer.get_data(
			args, FIG_3_EXPERIMENT_ID)
	make_fig3(data, environment_config)


if __name__ == '__main__':
	main()
