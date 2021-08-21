#! /usr/bin/env python
"""
Tools and analysis to debug charging problems.

TODO:
	- add additional sim options (aa_supply_in_charging, mechanistic_translation_supply) to match sim
"""

from __future__ import annotations

import argparse
import csv
import os
import pickle
from typing import Dict, Optional, Tuple
import webbrowser

import dash
import dash_core_components as dcc
import dash_html_components as html
import numpy as np
import plotly.graph_objs as go
import plotly.subplots

from models.ecoli.processes.metabolism import CONC_UNITS as METABOLISM_CONC_UNITS
from models.ecoli.processes.polypeptide_elongation import (calculate_trna_charging,
	CONC_UNITS, get_charging_params, get_ppgpp_params, ppgpp_metabolite_changes)
from wholecell.io.tablereader import TableReader
from wholecell.utils import constants, scriptBase, parallelization


# Set of charging/ppGpp parameters that should not be modified by a slider
CONSTANT_PARAMS = {'charging_mask', 'ppgpp_reaction_stoich', 'synthesis_index', 'degradation_index'}

PORT = 8050

# Element IDs
GRAPH_ID = 'graph-id'
LOW_TIMESTEP = 'low-timestep'
HIGH_TIMESTEP = 'high-timestep'

# Grid search columns
MEAN_COL = 'elongation rate mean (AA/s)'
STD_COL = 'elongation rate std'


def run_grid_search(charging, levels, index, search_params, timesteps):
	n_levels = len(levels)
	params = {}
	for i, (param, value) in enumerate(search_params.items()):
		params[param] = levels[index // n_levels**i % n_levels]

	v_ribs = []
	for timestep in timesteps:
		_, _, v, _, _ = charging.solve_timestep(timestep, param_adjustments=params)
		v_ribs.append(v)

	return v_ribs, params

class ChargingDebug(scriptBase.ScriptBase):
	def define_parameters(self, parser):
		super().define_parameters(parser)

		# Path args
		self.define_parameter_sim_dir(parser)
		self.define_parameter_variant_index(parser)
		parser.add_argument('-s', '--seed', type=int, default=0,
			help='The initial simulation number (int). The value will get'
				 ' formatted as a subdirectory name like "000000". Default = 0.')
		parser.add_argument('-g', '--generation', type=int, default=0,
			help='The generation number (int). The value will get formatted'
				 ' as a subdirectory name like "generation_000000". Default = 0.')
		parser.add_argument('-d', '--daughter', type=int, default=0,
			help='The daughter number (int). The value will get formatted as'
				 ' a subdirectory name like "000000". Default = 0.')

		# Sim options
		self.define_parameter_bool(parser, 'variable_elongation_translation',
			default_key='variable_elongation_translation',
			help='set if sims were run with variable_elongation_translation')

		# Debug options
		parser.add_argument('-c', '--cpus', type=int, default=1,
			help='Number of CPUs to use to run in parallel.')
		parser.add_argument('-o', '--output', default='charging-debug',
			help='Base name of output file to save.')
		parser.add_argument('--validation', type=int, default=1,
			help='Number of time steps to run for validation. If < 0, will run all.')
		parser.add_argument('--grid-search', action='store_true',
			help='If set, runs a grid search across parameters.')
		parser.add_argument('--grid-compare', nargs=2,
			help='Grid search output files to compare.')
		parser.add_argument('--interactive', action='store_true',
			help='If set, runs interactive analysis plots for debugging.')
		parser.add_argument('-p', '--port', type=int, default=PORT,
			help='The localhost port to use for the interactive webpage.')

	def update_args(self, args):
		super().update_args(args)

		# Extract data from args
		variant_dir_name, _, _ = args.variant_dir
		seed_str = '%06d' % (args.seed,)
		gen_str = 'generation_%06d' % (args.generation,)
		daughter_str = '%06d' % (args.daughter,)
		dirs = os.path.join(seed_str, gen_str, daughter_str)
		input_variant_directory = os.path.join(args.sim_path, variant_dir_name)

		# Set paths from args
		args.sim_data_file = os.path.join(input_variant_directory, 'kb',
			constants.SERIALIZED_SIM_DATA_MODIFIED)
		args.sim_out_dir = os.path.join(input_variant_directory, dirs, 'simOut')

	def load_data(self, sim_data_file: str, sim_out_dir: str) -> None:
		"""
		Loads sim_data and simulation output data from files and saves it as
		instance variables.

		Args:
			sim_data_file: path to the sim_data file for the simulation
			sim_out_dir: path to the simOut dir for the simulation
		"""

		with open(sim_data_file, 'rb') as f:
			self.sim_data = pickle.load(f)
		self.aa_from_trna = self.sim_data.process.transcription.aa_from_trna

		# Get charging parameters and separate to ones that will be adjusted or not
		charging_params = get_charging_params(self.sim_data,
			variable_elongation=self.variable_elongation)
		self.adjustable_charging_params = {
			param: value
			for param, value in charging_params.items()
			if param not in CONSTANT_PARAMS
			}
		self.constant_charging_params = {
			param: value
			for param, value in charging_params.items()
			if param in CONSTANT_PARAMS
			}

		# Get ppGpp parameters and separate to ones that will be adjusted or not
		ppgpp_params = get_ppgpp_params(self.sim_data)
		self.adjustable_ppgpp_params = {
			param: value
			for param, value in ppgpp_params.items()
			if param not in CONSTANT_PARAMS
			}
		self.constant_ppgpp_params = {
			param: value
			for param, value in ppgpp_params.items()
			if param in CONSTANT_PARAMS
			}

		# Listeners used
		growth_reader = TableReader(os.path.join(sim_out_dir, 'GrowthLimits'))
		kinetics_reader = TableReader(os.path.join(sim_out_dir, 'EnzymeKinetics'))
		main_reader = TableReader(os.path.join(sim_out_dir, 'Main'))

		# Load data
		self.synthetase_conc = CONC_UNITS * growth_reader.readColumn('synthetase_conc')[1:, :]
		self.uncharged_trna_conc = CONC_UNITS * growth_reader.readColumn('uncharged_trna_conc')[1:, :]
		self.charged_trna_conc = CONC_UNITS * growth_reader.readColumn('charged_trna_conc')[1:, :]
		self.aa_conc = CONC_UNITS * growth_reader.readColumn('aa_conc')[1:, :]
		self.ribosome_conc = CONC_UNITS * growth_reader.readColumn('ribosome_conc')[1:]
		self.fraction_aa_to_elongate = growth_reader.readColumn('fraction_aa_to_elongate')[1:, :]
		self.fraction_charged = growth_reader.readColumn('fraction_trna_charged')[1:, :]
		self.ppgpp_conc = CONC_UNITS * growth_reader.readColumn('ppgpp_conc')[1:]
		self.rela_conc = CONC_UNITS * growth_reader.readColumn('rela_conc')[1:]
		self.spot_conc = CONC_UNITS * growth_reader.readColumn('spot_conc')[1:]

		self.counts_to_molar = METABOLISM_CONC_UNITS * kinetics_reader.readColumn('countsToMolar')[1:]

		self.time_step_sizes = main_reader.readColumn('timeStepSec')[1:]
		self.n_time_steps = len(self.time_step_sizes)

	def solve_timestep(self,
			timestep: int,
			synthetase_adjustments: Optional[np.ndarray] = None,
			trna_adjustments: Optional[np.ndarray] = None,
			aa_adjustments: Optional[np.ndarray] = None,
			param_adjustments: Optional[Dict] = None,
			ribosome_adjustment: float = 1.,
			timestep_adjustment: float = 1.,
			) -> Tuple[np.ndarray, np.ndarray, float, int, int]:
		"""
		Calculates charging and elongation rate for a given timestep.

		Args:
			timestep: simulation timestep to select data from
			synthetase_adjustments: adjustments to scale synthetase concentrations
				up or down
			trna_adjustments: adjustments to scale tRNA concentrations
				up or down
			aa_adjustments: adjustments to scale amino acid concentrations
				up or down
			param_adjustments: adjustments to charging parameters
			ribosome_adjustment: adjustment to scale ribosome concentrations
			timestep_adjustment: adjustment to scale timesteps

		Returns:
			fraction_charged_per_trna: fraction charged for each tRNA
			fraction_charged: fraction charged of all tRNAs for each amino acid
			adjusted_v_rib: ribosome elongation rate (AA/s)
			n_synthesis: number of ppGpp synthesis reactions
			n_degradation: number of ppGpp degradation reactions

		TODO:
			Include adjustments for ppGpp, RelA, and SpoT concentrations
		"""

		n_aas = len(self.sim_data.molecule_groups.amino_acids)

		if synthetase_adjustments is None:
			synthetase_adjustments = np.ones(n_aas)
		if trna_adjustments is None:
			trna_adjustments = np.ones(n_aas)
		if aa_adjustments is None:
			aa_adjustments = np.ones(n_aas)
		if param_adjustments is None:
			param_adjustments = {}

		# Update parameters that need to be adjusted
		charging_params = {}
		for param, value in self.adjustable_charging_params.items():
			charging_params[param] = value * param_adjustments.get(param, 1)
		charging_params.update(self.constant_charging_params)
		ppgpp_params = {}
		for param, value in self.adjustable_ppgpp_params.items():
			ppgpp_params[param] = value * param_adjustments.get(param, 1)
		ppgpp_params.update(self.constant_ppgpp_params)

		# Adjust concentrations that are used in multiple locations
		ribosome_conc = self.ribosome_conc[timestep] * ribosome_adjustment
		uncharged_trna_conc = self.uncharged_trna_conc[timestep, :] * trna_adjustments
		charged_trna_conc = self.charged_trna_conc[timestep, :] * trna_adjustments
		f = self.fraction_aa_to_elongate[timestep, :]
		timestep_size = self.time_step_sizes[timestep] * timestep_adjustment

		# Calculate tRNA charging and resulting values
		fraction_charged, v_rib, _ = calculate_trna_charging(
			self.synthetase_conc[timestep, :] * synthetase_adjustments,
			uncharged_trna_conc,
			charged_trna_conc,
			self.aa_conc[timestep, :] * aa_adjustments,
			ribosome_conc,
			f,
			charging_params,
			time_limit=timestep_size,
			)
		fraction_charged_per_aa = fraction_charged @ self.aa_from_trna
		adjusted_v_rib = v_rib / ribosome_conc.asNumber(CONC_UNITS)

		# Update tRNA concentrations to reflect charging
		total_trna_conc = uncharged_trna_conc + charged_trna_conc
		updated_uncharged_trna_conc = total_trna_conc * (1 - fraction_charged)
		updated_charged_trna_conc = total_trna_conc * fraction_charged

		# Calculate ppGpp reaction rates
		_, n_synthesis, n_degradation, _, _, _ = ppgpp_metabolite_changes(
			updated_uncharged_trna_conc,
			updated_charged_trna_conc,
			ribosome_conc,
			f,
			self.rela_conc[timestep],
			self.spot_conc[timestep],
			self.ppgpp_conc[timestep],
			self.counts_to_molar[timestep],
			v_rib,
			charging_params,
			ppgpp_params,
			timestep_size,
			)

		return fraction_charged_per_aa, fraction_charged, adjusted_v_rib, n_synthesis, n_degradation

	def validation(self, n_steps: int) -> None:
		"""
		Performs a validation check to makes sure solving the model from
		loaded data matches the objective from the original solution during
		the simulation.

		Args:
			n_steps: number of timesteps to check
				if 0: does not check
				if <0: runs all timepoints from the simulation
		"""

		if n_steps == 0:
			return
		elif n_steps < 0:
			n_steps = self.n_time_steps
		else:
			n_steps = min(self.n_time_steps, n_steps)

		print('Running validation to check output...')
		for timestep in range(n_steps):
			fraction_charged, _, _, _, _ = self.solve_timestep(timestep)
			if np.any(fraction_charged != self.fraction_charged[timestep]):
				raise ValueError(f'Charging fraction does not match for time step {timestep}')
		print('All {} timesteps match the results from the whole-cell model.'.format(n_steps))

	def grid_search(self, filename, n_steps=20, n_levels=9, low_level=-1, high_level=1, cpus=1):
		"""
		Perform a grid search on the charging parameters to determine the effect
		on elongation rate.  Results are saved to a tsv file that can be used
		with plot_grid_results to compare two different grid search results.

		TODO:
			accept args from the command line
		"""

		filename = f'{filename}.tsv'

		levels = np.logspace(low_level, high_level, n_levels)
		search_params = {k: v for k, v in self.adjustable_charging_params.items() if k != 'max_elong_rate'}
		n_params = len(search_params)
		timesteps = np.arange(self.n_time_steps)[::int(np.ceil(self.n_time_steps / n_steps))]
		n_samples = n_levels**n_params

		with open(filename, 'w') as f:
			writer = csv.DictWriter(f, fieldnames=list(search_params.keys()) + [MEAN_COL, STD_COL], delimiter='\t')
			writer.writeheader()

			def callback(result):
				v_ribs, params = result
				writer.writerow({MEAN_COL: np.mean(v_ribs), STD_COL: np.std(v_ribs), **params})

			# Run timesteps in parallel
			pool = parallelization.pool(num_processes=cpus)
			results = [
				pool.apply_async(run_grid_search, (self, levels, index, search_params, timesteps), callback=callback)
				for index in range(n_samples)
				]
			pool.close()
			pool.join()

			# Check for errors
			for result in results:
				if not result.successful():
					result.get()

	def plot_grid_results(self, filename: str, path1: str, path2: str):
		"""
		Compares results from two different grid searches in an interactive
		plot.  Useful for checking sets of parameters that give desired results
		in different simulations.

		Args:
			filename: path to the html file with the comparison plot
			path1: path to the tsv file results from one grid search
			path2: path to the tsv file results from another grid search
		"""

		def load_data(path):
			data = {}
			with open(path) as f:
				reader = csv.reader(f, delimiter='\t')
				headers = next(reader)
				mean_col = headers.index(MEAN_COL)
				std_col = headers.index(STD_COL)
				for line in reader:
					data[tuple(zip(headers, line[:mean_col]))] = {
						'mean': float(line[mean_col]),
						'std': float(line[std_col]),
						}

			return data

		filename = f'{filename}.html'

		data1 = load_data(path1)
		data2 = load_data(path2)

		labels = list(data1.keys() | data2.keys())
		x = np.array([data1.get(label, {}).get('mean', 0) for label in labels])
		y = np.array([data2.get(label, {}).get('mean', 0) for label in labels])
		error_x = dict(
			type='data',
			array=[data1.get(label, {}).get('std', 0) for label in labels],
			visible=False,  # Need to find better way of displaying before making visible
			)
		error_y = dict(
			type='data',
			array=[data2.get(label, {}).get('std', 0) for label in labels],
			visible=False,  # Need to find better way of displaying before making visible
			)

		fig = go.Figure(data=go.Scatter(x=x, y=y, error_x=error_x, error_y=error_y, mode='markers', text=labels))
		fig.write_html(filename)

	def interactive_debug(self, port: int):
		"""
		Run an interactive app in a browser to debug charging.
		"""

		app = self.create_app()
		webbrowser.open_new(f'http://127.0.0.1:{port}/')
		app.run_server(port=port)

	def create_app(self) -> dash.Dash:
		"""
		Create the Dash app to run in the browser for interactive mode.
		"""

		app = dash.Dash()

		charging_param_ids = sorted(self.adjustable_charging_params)
		ppgpp_param_ids = sorted(self.adjustable_ppgpp_params)
		all_param_ids = charging_param_ids + ppgpp_param_ids
		aa_ids = self.sim_data.molecule_groups.amino_acids
		n_aas = len(aa_ids)

		# Slider elements
		def add_slider(id_):
			text_id = f'{id_}-text'
			slider_id = f'{id_}-slider'
			div = html.Div(style={'display': 'grid', 'grid-template-columns': '10% 90%'}, children=[
				html.Plaintext('1.000', id=text_id),
				dcc.Slider(id=slider_id, **slider_options),
				])
			@app.callback(dash.dependencies.Output(text_id, 'children'),
				[dash.dependencies.Input(slider_id, 'value')])
			def display_value(value):
				return f'{10**value:.3f}'

			return div

		slider_options = dict(value=0, min=-2, max=2, step=0.01, marks={i: {'label': 10**i} for i in range(-2, 3)})
		aa_slider_style = {'display': 'grid', 'grid-template-columns': '15% 25% 25% 25%'}
		aa_slider_headers = [html.Div(style=aa_slider_style, children=[
			html.Plaintext(''),
			html.Plaintext('Synthetase concentrations', style={'text-align': 'center'}),
			html.Plaintext('tRNA concentrations', style={'text-align': 'center'}),
			html.Plaintext('Amino acid concentrations', style={'text-align': 'center'}),
			])]
		aa_sliders = [
			html.Div(style=aa_slider_style, children=[
				html.Plaintext(f'{aa[:-3]}:'),
				add_slider(f'{aa}-synthetase'),
				add_slider(f'{aa}-trna'),
				add_slider(f'{aa}-aa'),
				])
			for aa in aa_ids
			]
		other_sliders_style = {'display': 'grid', 'grid-template-columns': '30% 70%'}
		charging_param_headers = [html.Div(style=other_sliders_style, children=[
			html.Plaintext(''),
			html.Plaintext('Charging parameters', style={'text-align': 'center'}),
			])]
		charging_param_sliders = [
			html.Div(style=other_sliders_style, children=[
				html.Plaintext(f'{param}:'),
				add_slider(f'{param}-slider'),
				])
			for param in charging_param_ids
			]
		ppgpp_param_headers = [html.Div(style=other_sliders_style, children=[
			html.Plaintext(''),
			html.Plaintext('ppGpp parameters', style={'text-align': 'center'}),
			])]
		ppgpp_param_sliders = [
			html.Div(style=other_sliders_style, children=[
				html.Plaintext(f'{param}:'),
				add_slider(f'{param}-slider'),
				])
			for param in ppgpp_param_ids
			]
		other_headers = [html.Div(style=other_sliders_style, children=[
			html.Plaintext(''),
			html.Plaintext('Other inputs', style={'text-align': 'center'}),
			])]
		other_sliders = [
			html.Div(style=other_sliders_style, children=[
				html.Plaintext('Ribosome conc:'),
				add_slider('ribosome-slider'),
				]),
			html.Div(style=other_sliders_style, children=[
				html.Plaintext('Timestep:'),
				add_slider('timestep-slider'),
				]),
			]

		sliders = html.Div(
			style={'maxHeight': '500px', 'overflow': 'scroll', 'display': 'grid', 'grid-template-columns': '70% 30%'},
			children=[
				html.Div(children=aa_slider_headers + aa_sliders),
				html.Div(children=charging_param_headers
					+ charging_param_sliders
					+ ppgpp_param_headers
					+ ppgpp_param_sliders
					+ other_headers
					+ other_sliders),
			])

		# Slider inputs
		synthetase_inputs = [
			dash.dependencies.Input(f'{aa}-synthetase-text', 'children')
			for aa in aa_ids
			]
		trna_inputs = [
			dash.dependencies.Input(f'{aa}-trna-text', 'children')
			for aa in aa_ids
			]
		aa_inputs = [
			dash.dependencies.Input(f'{aa}-aa-text', 'children')
			for aa in aa_ids
			]
		param_inputs = [
			dash.dependencies.Input(f'{param}-slider-text', 'children')
			for param in all_param_ids
			]

		# Page layout
		app.layout = html.Div(children=[
			html.H2('Charging debugger'),
			html.Plaintext('Timestep limits (lower, upper):'),
			dcc.Input(id=LOW_TIMESTEP, type='number', value=0),
			dcc.Input(id=HIGH_TIMESTEP, type='number', value=10),
			sliders,
			dcc.Graph(id=GRAPH_ID, style={'height': '750px'})
			])

		# Register callback to update plot when selections change
		# First arg for Output/Input selects the page object
		# Second arg for Output/Input sets or gets a kwarg from the dcc function
		@app.callback(
			dash.dependencies.Output(GRAPH_ID, 'figure'),
			[
				dash.dependencies.Input(LOW_TIMESTEP, 'value'),
				dash.dependencies.Input(HIGH_TIMESTEP, 'value'),
				dash.dependencies.Input('ribosome-slider-text', 'children'),
				dash.dependencies.Input('timestep-slider-text', 'children'),
				*synthetase_inputs,
				*trna_inputs,
				*aa_inputs,
				*param_inputs,
			])
		def update_graph(
				init_t: int, final_t: int,
				ribosome_adjustment: str, timestep_adjustment: str,
				*param_inputs: str,
				) -> Dict:
			"""
			Update the plot based on selection changes.

			Returns:
				plotly figure
			"""

			ribosome_adjustment = float(ribosome_adjustment)
			timestep_adjustment = float(timestep_adjustment)
			synthetase_adjustments = np.array(param_inputs[:n_aas], float)
			trna_adjustments = np.array(param_inputs[n_aas:2*n_aas], float)
			aa_adjustments = np.array(param_inputs[2*n_aas:3*n_aas], float)
			param_adjustments = {param: float(value) for param, value in zip(all_param_ids, param_inputs[3*n_aas:])}

			v_rib = []
			f_charged = []
			n_synth = []
			n_deg = []
			t = np.arange(init_t, final_t)
			for timestep in t:
				_, f, v, synth, deg = self.solve_timestep(
					timestep,
					synthetase_adjustments=synthetase_adjustments,
					trna_adjustments=trna_adjustments,
					aa_adjustments=aa_adjustments,
					param_adjustments=param_adjustments,
					ribosome_adjustment=ribosome_adjustment,
					timestep_adjustment=timestep_adjustment,
					)
				v_rib.append(v)
				f_charged.append(f)
				n_synth.append(synth)
				n_deg.append(deg)

			f_charged = np.array(f_charged).T
			net_ppgpp = np.array(n_synth) - np.array(n_deg)

			fig = plotly.subplots.make_subplots(rows=1, cols=3)
			fig.append_trace(go.Scatter(x=t, y=v_rib, name='Elongation rate'), row=1, col=1)
			for f, aa in zip(f_charged, aa_ids):
				fig.append_trace(go.Scatter(x=t, y=f, name=aa), row=1, col=2)
			fig.append_trace(go.Scatter(x=t, y=n_synth, name='ppGpp synthesis reactions'), row=1, col=3)
			fig.append_trace(go.Scatter(x=t, y=n_deg, name='ppGpp degradation reactions'), row=1, col=3)
			fig.append_trace(go.Scatter(x=t, y=net_ppgpp, name='Net ppGpp reactions'), row=1, col=3)

			fig.update_xaxes(title_text='Timestep', row=1, col=2)
			fig.update_yaxes(title_text='Elongation rate (AA/s)', row=1, col=1)
			fig.update_yaxes(title_text='Fraction charged', row=1, col=2)
			fig.update_yaxes(title_text='Number of ppGpp reaction', row=1, col=3)

			return fig

		return app

	def run(self, args: argparse.Namespace) -> None:
		self.variable_elongation = args.variable_elongation_translation
		self.load_data(args.sim_data_file, args.sim_out_dir)
		self.validation(args.validation)
		if args.grid_search:
			self.grid_search(args.output, cpus=args.cpus)
		if args.grid_compare:
			self.plot_grid_results(args.output, *args.grid_compare)
		if args.interactive:
			self.interactive_debug(args.port)


if __name__ == '__main__':
	analysis = ChargingDebug()
	analysis.cli()
