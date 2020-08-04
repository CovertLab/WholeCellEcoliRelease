'''Generate figures for wcEcoli colony simulation

For usage information, run:
	python make_figures.py -h
'''
import argparse
import os

from vivarium.analysis.analyze import Analyzer
from vivarium.plots.multibody_physics import plot_tags


TAG_PATH_NAME_MAP = {
	('boundary', 'bulk_molecules_report', 'XAPB-MONOMER[i]'): 'XapB',
	(
		'boundary', 'bulk_molecules_report', 'TRANS-CPLX-201[s]'
	): 'AcrAB-TolC',
}
OUT_DIR = os.path.join('environment', 'figs')
FILE_EXTENSION = 'eps'
FIG_2_EXPERIMENT_ID = '20200731.144039'


def make_fig2(data, environment_config):
	tags_data = Analyzer.format_data_for_tags(data, environment_config)
	plot_config = {
		'out_dir': OUT_DIR,
		'tagged_molecules': TAG_PATH_NAME_MAP.keys(),
		'filename': 'fig2.{}'.format(FILE_EXTENSION),
		'tag_path_name_map': TAG_PATH_NAME_MAP,
		'tag_label_size': 36,
	}
	plot_tags(tags_data, plot_config)


def main():
	if not os.path.exists(OUT_DIR):
		os.makedirs(OUT_DIR)
	parser = argparse.ArgumentParser()
	Analyzer.add_connection_args(parser)
	args = parser.parse_args()
	data, environment_config = Analyzer.get_data(
		args, FIG_2_EXPERIMENT_ID)
	make_fig2(data, environment_config)


if __name__ == '__main__':
    main()
