'''Generate figures for wcEcoli colony simulation

For usage information, run:
	python make_figures.py -h
'''
import argparse
import os

from vivarium.analysis.analyze import Analyzer
from vivarium.plots.multibody_physics import plot_tags, plot_snapshots


TAG_PATH_NAME_MAP = {
	('boundary', 'bulk_molecules_report', 'XAPB-MONOMER[i]'): 'XapB',
	(
		'boundary', 'bulk_molecules_report', 'TRANS-CPLX-201[s]'
	): 'AcrAB-TolC',
}
OUT_DIR = os.path.join('environment', 'figs')
FILE_EXTENSION = 'eps'
FIG_2_EXPERIMENT_ID = '20200805.012237'
FIG_3_EXPERIMENT_ID = '20200805.012237'


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
	parser = argparse.ArgumentParser()
	Analyzer.add_connection_args(parser)
	args = parser.parse_args()

	data_2, environment_config_2 = Analyzer.get_data(
		args, FIG_2_EXPERIMENT_ID)
	make_fig2(data_2, environment_config_2)

	data_3, environment_config_3 = Analyzer.get_data(
		args, FIG_3_EXPERIMENT_ID)
	make_fig3(data_3, environment_config_3)


if __name__ == '__main__':
	main()
