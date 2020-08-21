'''Generate figures for wcEcoli colony simulation

For usage information, run:
	python make_figures.py -h
'''
import argparse
import os
import sys

from vivarium.analysis.analyze import Analyzer
from vivarium.plots.colonies import plot_metric_across_experiments
from vivarium.plots.multibody_physics import plot_tags, plot_snapshots
from vivarium.analysis.expression_survival_dotplot import (
	plot_expression_survival)
from vivarium.core.composition import plot_agents_multigen

import wholecell.utils.filepath as fp


PUMP_PATH = (
	'boundary', 'bulk_molecule_concentrations', 'TRANS-CPLX-201[s]')
TAG_PATH_NAME_MAP = {
	('boundary', 'bulk_molecules_report', 'XAPB-MONOMER[i]'): 'XapB',
	(
		'boundary', 'bulk_molecules_report', 'TRANS-CPLX-201[s]'
	): 'AcrAB-TolC',
}
COLONY_MASS_PATH = ('mass',)
OUT_DIR = os.path.join('environment', 'figs')
FILE_EXTENSION = 'pdf'
FIG_2_EXPERIMENT_ID = '20200820.202016'
FIG_3_EXPERIMENT_ID = FIG_2_EXPERIMENT_ID
FIG_4A_EXPERIMENT_ID = FIG_2_EXPERIMENT_ID
FIG_4B_EXPERIMENT_ID = '20200820.235622'
FIG_4_5_EXPERIMENT_IDS = {
	'0.002 mM': '20200818.174841',
	'0.011 mM': '20200818.215147',
	'0.0155 mM': '20200819.151444',
	'0.01775 mM': '20200819.175108',
	'0.018875 mM': '20200819.203802',
	'0.02 mM': '20200817.224609',
}
FIG_6_EXPERIMENT_ID = '20200817.224609'
FIG_7_EXPERIMENT_ID = FIG_6_EXPERIMENT_ID
FIG_6_TIME_RANGE = (0.5, 1)
METADATA_FILE = 'metadata.json'


def get_metadata():
	metadata = {
		'git_hash': fp.run_cmdline('git rev-parse HEAD'),
		'git_branch': fp.run_cmdline('git symbolic-ref --short HEAD'),
		'time': fp.timestamp(),
		'python': sys.version.splitlines()[0],
		'fig2_id': FIG_2_EXPERIMENT_ID,
		'fig3_id': FIG_3_EXPERIMENT_ID,
		'fig4a_id': FIG_4A_EXPERIMENT_ID,
		'fig4b_id': FIG_4B_EXPERIMENT_ID,
		'fig6': FIG_6_EXPERIMENT_ID,
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
		'tag_label_size': 48,
	}
	plot_tags(tags_data, plot_config)


def make_snapshots_figure(data, environment_config, name):
	snapshots_data = Analyzer.format_data_for_snapshots(
		data, environment_config)
	plot_config = {
		'out_dir': OUT_DIR,
		'filename': '{}.{}'.format(name, FILE_EXTENSION),
		'include_fields': ['nitrocefin'],
		'field_label_size': 48,
	}
	plot_snapshots(snapshots_data, plot_config)


def make_fig4_5(data_dict):
	path_ts_dict = {
		key: Analyzer.format_data_for_colony_metrics(value)
		for key, value in data_dict.items()
	}
	fig = plot_metric_across_experiments(
		path_ts_dict, COLONY_MASS_PATH, ylabel='Colony Mass (mg)')
	fig.savefig(os.path.join(
		OUT_DIR, 'fig4_5.{}'.format(FILE_EXTENSION)))


def make_fig6(data):
	fig = plot_expression_survival(
		data, PUMP_PATH,
		'Average AcrAB-TolC Concentration (mmol/L) Over Cell Lifetime',
		FIG_6_TIME_RANGE,
	)
	fig.savefig(os.path.join(
		OUT_DIR, 'fig6.{}'.format(FILE_EXTENSION)))


def make_fig7(data):
	settings = {
		'include_paths': [PUMP_PATH],
		'titles_map': {
			PUMP_PATH: 'AcrAB-TolC Concentration',
		},
		'ylabels_map': {
			PUMP_PATH: 'Concentration (mM)',
		},
	}
	plot_agents_multigen(
		data, settings, OUT_DIR, 'fig7.{}'.format(FILE_EXTENSION))


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
	make_snapshots_figure(data, environment_config, 'fig3')

	if FIG_4A_EXPERIMENT_ID != FIG_3_EXPERIMENT_ID:
		data, environment_config = Analyzer.get_data(
			args, FIG_4A_EXPERIMENT_ID)
	make_snapshots_figure(data, environment_config, 'fig4a')

	data, environment_config = Analyzer.get_data(
		args, FIG_4B_EXPERIMENT_ID)
	make_snapshots_figure(data, environment_config, 'fig4b')

	data_dict = dict()
	for key, exp_id in FIG_4_5_EXPERIMENT_IDS.items():
		exp_data, _ = Analyzer.get_data(args, exp_id)
		data_dict[key] = exp_data
	make_fig4_5(data_dict)

	data, _ = Analyzer.get_data(args, FIG_6_EXPERIMENT_ID)
	make_fig6(data)

	if FIG_6_EXPERIMENT_ID != FIG_7_EXPERIMENT_ID:
		data, _ = Analyzer.get_data(args, FIG_7_EXPERIMENT_ID)
	make_fig7(data)


if __name__ == '__main__':
	main()
