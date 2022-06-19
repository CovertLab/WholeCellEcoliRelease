#! /usr/bin/env python

"""
Explore simulation data in an interactive manner with the ability to select
datasets and graph types.

TODO:
	- default column option? - start with time as x (select from command line)
	- add reducing options - downsampling
	- access sim_data/validation_data arrays
	- select multiple y datasets
	- add multigen/cohort/variant selection options
	- add '*' selection option (select all) when multiple directories are found
"""

from __future__ import annotations

import argparse
import os
import re
from typing import Dict, List, Optional, Set, Tuple, Union
import webbrowser

import dash
import dash_core_components as dcc
import dash_html_components as html
import numpy as np
import plotly.graph_objs as go

from wholecell.io.tablereader import DoesNotExistError, TableReader
from wholecell.utils import scriptBase
import wholecell.utils.filepath as fp


PORT = 8050

# Options from drop down menu
PLOT_OPTIONS = {
	'line': {'function': go.Scatter},
	'bar': {'function': go.Bar, 'layout_options': {'barmode': 'stack'}},  # this could use more testing
	}

# Object IDs
GRAPH_ID = 'graph'
PLOT_SELECTION = 'plot-selector'
DATA_SELECTION_ID = 'Data selection:'
X_DATA_OPTIONS_ID = 'x-data-options'
Y_DATA_OPTIONS_ID = 'y-data-options'
ADD_X_ID = 'Update x'
ADD_Y_ID = 'Update y'
BUTTON_VALUE_TEMPLATE = '{} value'
SEPARATOR = '<>'
VALUE_JOIN = f'{{}}{SEPARATOR}{{}}'
DATA_OPTIONS = ['mean', 'normalized', 'log']


def get_vals(d: Dict, k: Union[str, List[str]]):
	"""Get values from a nested dictionary with nested keys provided."""

	if isinstance(k, str):
		k = k.split(SEPARATOR)

	if len(k) > 1:
		return get_vals(d[k[0]], k[1:])
	else:
		return d.get(k[0])

def load_listener(selection: str) -> Tuple[np.ndarray, List[str]]:
	"""
	Load data from a listener specified by the drop down selection directory
	string.

	Args:
		selection: directory path to simOut with listener and column as last two
			values (all separated by SEPARATOR)

	Returns:
		data: 2D data array (n timepoints, m data series)
		labels: labels corresponding to each m data series
	"""

	split = selection.split(SEPARATOR)
	path = os.path.join(*split[:-2], 'simOut')
	listener = split[-2]
	column = split[-1]

	reader = TableReader(os.path.join(path, listener))
	data = reader.readColumn(column, squeeze=False)

	# Read subcolumn attributes
	try:
		subcolumns = reader.readAttribute('subcolumns')
	except DoesNotExistError:
		subcolumns = {}
	if column in subcolumns:
		labels = reader.readAttribute(subcolumns[column])
	else:
		labels = list(range(data.shape[1]))

	return data, labels

def create_app(data_structure: Dict) -> dash.Dash:
	"""
	Create the dash app to serve the webpage and content.

	Args:
		data_structure: nested directory structure of possible simulations
			(see parse_data_structure())

	Returns:
		app: dash app that can run a server for the interactive plot
	"""

	def data_selection(
			app: dash.Dash,
			data_structure: Dict,
			id_: str,
			defaults: Optional[Set[str]] = None,
			multi: bool = False,
			) -> Tuple[html.Div, dash.dependencies.State]:
		"""
		Create div to hold drop down menus for data selection.

		Args:
			app: dash app
			data_structure: nested directory structure of possible simulations
				(see parse_data_structure())
			id_: ID for the drop down object on the webpage
			defaults: default values to start with in a drop down if any
				selection options match a value in this set
			multi: if True, creates drop down menus with the option of making
				multiple selections

		Returns:
			div: div containing all drop down menus
			value: value for the bottom drop down menu for the path to data
				selected

		TODO:
			- get multi selection working
			- handle initialization for x and y separately
		"""

		def get_selection_options(
				data_structure: Dict,
				parent_value: str,
				defaults: Set[str],
				current: Optional[str] = None,
				) -> Tuple[List[Dict[str, str]], Optional[str]]:
			"""
			Get the selection options for drop down menus based on the parent
			drop down value indexed into the data structure.

			Args:
				data_structure: nested directory structure of possible
					simulations (see parse_data_structure())
				parent_value: value of the parent drop down menu
				defaults: default values to start with in a drop down if any
					selection options match a value in this set
				current: current value of the drop down, if set, keeps this
					value if it is still a valid selection option

			Returns:
				options: drop down menu options with a display label and
					correponding value
				value: the default value to select for the drop down menu
			"""

			# Empty return values in case path or vals are not specified
			options = []
			value = None

			if parent_value is not None:
				vals = get_vals(data_structure, parent_value)
				if vals is not None:
					for val in sorted(vals):
						joined = VALUE_JOIN.format(parent_value, val)
						options.append({
							'label': val,
							'value': joined,
							})
						if val in defaults or value is None:
							value = joined

					if current is not None:
						current = current.split(SEPARATOR)[-1]
						if current in vals:
							value = VALUE_JOIN.format(parent_value, current)

			return options, value

		def add_children(
				children: List,
				base_id: str,
				parent_id: str,
				parent_value: str,
				data_structure: Dict,
				defaults: Set[str],
				count: int = 0
				) -> Tuple[List, int]:
			"""
			Recursively add drop down selections dependent on the parent
			drop down menu.  Registers a callback to update the drop down when
			the parent is updated.

			Args:
				children: children for the HTML div (header and drop down menus)
				base_id: ID for the top level drop down object on the webpage
				parent_id: ID for the parent drop down object on the webpage
				parent_value: value for the parent drop down
				data_structure: nested directory structure of possible simulations
					(see parse_data_structure())
				defaults: default values to start with in a drop down if any
					selection options match a value in this set
				count: number of children drop down menus added

			Returns:
				children: children for the HTML div (header and drop down menus)
					with new drop down menus added
				count: number of new drop downs added
			"""

			options, value = get_selection_options(data_structure, parent_value, defaults)

			if value is None:
				# Bottom of the data structure - no more selection is needed
				return children, count
			else:
				# New drop down menu for nested selection
				count += 1
				sub_id = f'{base_id}{count}'
				children.append(dcc.Dropdown(
					id=sub_id,
					multi=multi,
					options=options,
					value=value,
					))

				# Register callback to update list options when parent changes
				@app.callback(
					[
						dash.dependencies.Output(sub_id, 'options'),
						dash.dependencies.Output(sub_id, 'value'),
					],
					[dash.dependencies.Input(parent_id, 'value')],
					[dash.dependencies.State(sub_id, 'value')],
					prevent_initial_call=True)
				def update(parent_value: str, current: str) -> Tuple[List[Dict[str, str]], Optional[str]]:
					"""Update valid selection based on the parent value"""
					return get_selection_options(data_structure, parent_value, defaults, current=current)

				return add_children(children, base_id, sub_id, value, data_structure, defaults, count=count)

		if defaults is None:
			defaults = set()

		value = next(iter(data_structure))

		# Create children for the drop down div
		children = [
			html.H2(id_),
			dcc.Dropdown(
				id=id_,
				options=[{
					'label': os.path.basename(d),
					'value': d,
					} for d in data_structure],
				value=value,
				)
			]
		children, n_added = add_children(children, id_, id_, value, data_structure, defaults)

		div = html.Div(children=children)
		input_value = dash.dependencies.State(f'{id_}{n_added}', 'value')

		return div, input_value

	# Create webpage layout
	app = dash.Dash()
	input_div, input_value = data_selection(app, data_structure, DATA_SELECTION_ID, defaults={'Main', 'time'})
	app.layout = html.Div(children=[
		html.H1('Whole-cell simulation explorer'),
		html.Div(children=[
			html.H2('Plot selection:'),
			dcc.Dropdown(
				id=PLOT_SELECTION,
				options=[{
					'label': o,  # display name
					'value': o,  # value passed through callback, must be str
					} for o in PLOT_OPTIONS],
				value=next(iter(PLOT_OPTIONS)),
				),
			input_div,
			html.Div(children=[
				html.Plaintext('x data options: '),
				dcc.Checklist(
					id=X_DATA_OPTIONS_ID,
					options=[{
						'label': o,  # display name
						'value': o,  # value passed through callback, must be str
						} for o in DATA_OPTIONS],
					value=[],
					),
				]),
			html.Div(children=[
				html.Plaintext('y data options: '),
				dcc.Checklist(
					id=Y_DATA_OPTIONS_ID,
					options=[{
						'label': o,  # display name
						'value': o,  # value passed through callback, must be str
						} for o in DATA_OPTIONS],
					value=[],
					),
				]),
			html.Div(children=[
				html.Button(ADD_X_ID, id=ADD_X_ID),
				html.Button(ADD_Y_ID, id=ADD_Y_ID),
				]),
			html.Div(children=[
				html.Plaintext('x: ', id=BUTTON_VALUE_TEMPLATE.format(ADD_X_ID)),
				html.Plaintext('y: ', id=BUTTON_VALUE_TEMPLATE.format(ADD_Y_ID)),
				]),
			]),
		dcc.Graph(id=GRAPH_ID, style={'height': '600px'}),
		])

	# Only update axis values on button click
	for button in [ADD_X_ID, ADD_Y_ID]:
		@app.callback(
			[
				dash.dependencies.Output(button, 'value'),
				dash.dependencies.Output(BUTTON_VALUE_TEMPLATE.format(button), 'children'),
			],
			dash.dependencies.Input(button, 'n_clicks'),  # needed for callback trigger
			[
				input_value,
				dash.dependencies.State(BUTTON_VALUE_TEMPLATE.format(button), 'children'),
			])
		def update_axis(n_clicks, val, previous):
			new_text = ' '.join(previous.split(' ')[:1] + val.split(SEPARATOR))
			return val, new_text

	# Register callback to update plot when selections change
	# First arg for Output/Input selects the page object
	# Second arg for Output/Input sets or gets a kwarg from the dcc function
	@app.callback(
		dash.dependencies.Output(GRAPH_ID, 'figure'),
		[
			dash.dependencies.Input(PLOT_SELECTION, 'value'),
			dash.dependencies.Input(ADD_X_ID, 'value'),
			dash.dependencies.Input(ADD_Y_ID, 'value'),
			dash.dependencies.Input(X_DATA_OPTIONS_ID, 'value'),
			dash.dependencies.Input(Y_DATA_OPTIONS_ID, 'value'),
		])
	def update_graph(plot_id: str, x_input: str, y_input: str, x_options: List[str], y_options: List[str]) -> Dict:
		"""
		Update the plot based on selection changes.

		Args:
			plot_id: ID for the plot type in PLOT_OPTIONS to display
			x_input: directory path to simOut with listener and column as last
				two values (all separated by SEPARATOR) for x data
			y_input: directory path to simOut with listener and column as last
				two values (all separated by SEPARATOR) for y data
			x_options: options from check boxes to apply transformations to the
				x data
			y_options: options from check boxes to apply transformations to the
				y data

		Returns:
			plotly figure dict
		"""

		def adjust_data(data, options):
			if data.shape[1] > 1:
				if 'normalized' in options:
					data /= data[1, :]  # use 1 as first index since a lot of listeners start at 0 for first entry
				if 'mean' in options:
					data = data.mean(0).reshape(1, -1)  # need to keep as 2D array
				if 'log' in options:
					data = np.log10(data)
			return data

		if x_input is None or y_input is None:
			return {}

		x_data, x_labels = load_listener(x_input)
		y_data, y_labels = load_listener(y_input)

		x_data = adjust_data(x_data, x_options)
		y_data = adjust_data(y_data, y_options)

		plot = PLOT_OPTIONS[plot_id]
		plot_options = plot.get('plot_options', {})
		layout_options = plot.get('layout_options', {})
		if x_data.shape == y_data.shape:
			traces = [
				plot['function'](x=x_data[:, idx], y=y_data[:, idx], name=col, **plot_options)
				for idx, col in enumerate(y_labels)
				]
		else:
			traces = [
				plot['function'](x=x_data[:, 0], y=y_data[:, idx], name=col, **plot_options)
				for idx, col in enumerate(y_labels)
				]

		# Dict used to update 'figure' for dcc.Graph object GRAPH_ID
		return {
			'data': traces,
			'layout': go.Layout(
				xaxis_title=x_input.split(SEPARATOR)[-1],
				yaxis_title=y_input.split(SEPARATOR)[-1],
				hovermode='closest',
				**layout_options),
			}

	return app

class AnalysisInteractive(scriptBase.ScriptBase):
	def define_parameters(self, parser):
		super().define_parameters(parser)

		self.define_parameter_sim_dir(parser, default=os.path.join(fp.ROOT_PATH, 'out'))

	def parse_data_structure(self, path: str) -> Dict:
		"""
		Parses the directory structure for simulation data that can be loaded.

		Args:
			path: path to a directory containing different simulation sets
				(eg. out/) or a path to a single simulation output directory
				(eg. out/single-simulation/)

		Returns:
			experiments: nested dictionary including paths to simOut and
				listeners and columns for all simulation directories found
		"""

		experiments = {}  # type: Dict

		# Look for variants in current directory or one directory deep
		variant_dirs = [v[0] for v in self.list_variant_dirs(path)]
		if len(variant_dirs):
			experiments[path] = {d: {} for d in variant_dirs}
		else:
			for directory in sorted(os.listdir(path)):
				sim_path = os.path.join(path, directory)
				if not os.path.isdir(sim_path):
					continue

				variant_dirs = [v[0] for v in self.list_variant_dirs(sim_path)]
				if len(variant_dirs):
					experiments[sim_path] = {d: {} for d in variant_dirs}

		# Find all possible simulations to select
		# TODO: more efficient or cleaner way of doing this
		# TODO: only find matches when selected on the page instead of prior
		found_listeners = False
		seed_regex = re.compile('[0-9]{6}')
		generation_regex = re.compile('generation_[0-9]{6}')
		daughter_regex = re.compile('[0-9]{6}')
		for base, variants in experiments.items():
			for variant, seed_dict in variants.items():
				variant_dir = os.path.join(base, variant)
				for seed in os.listdir(variant_dir):
					if seed_regex.match(seed):
						gen_dict = seed_dict.get(seed, {})
						seed_dir = os.path.join(variant_dir, seed)
						for gen in os.listdir(seed_dir):
							if generation_regex.match(gen):
								daughter_dict = gen_dict.get(gen, {})
								gen_dir = os.path.join(seed_dir, gen)
								for daughter in os.listdir(gen_dir):
									if daughter_regex.match(daughter):
										listener_dict = daughter_dict.get(daughter, {})
										sim_out_dir = os.path.join(gen_dir, daughter, 'simOut')
										for listener in os.listdir(sim_out_dir):
											column_dict = listener_dict.get(listener, {})
											listener_dir = os.path.join(sim_out_dir, listener)
											if os.path.isdir(listener_dir):
												for column in os.listdir(listener_dir):
													if '.' not in column:  # TODO: handle .cPickle columns
														column_dict[column] = None
														found_listeners = True
												listener_dict[listener] = column_dict
										daughter_dict[daughter] = listener_dict
								gen_dict[gen] = daughter_dict
						seed_dict[seed] = gen_dict

		if len(experiments) == 0 or not found_listeners:
			raise ValueError(f'Could not find valid simulations in "{path}"'
				' or its immediate subdirectories. Make sure the provided path'
				' is to a top level simulation output directory or directory'
				' containing one or more simulation output directories.')

		return experiments

	def run(self, args: argparse.Namespace) -> None:
		data_structure = self.parse_data_structure(args.sim_dir)
		app = create_app(data_structure)

		# Serve interactive page (may take a second to load, reload if necessary)
		webbrowser.open_new(f'http://127.0.0.1:{PORT}/')
		app.run_server(port=PORT)


if __name__ == '__main__':
	analysis = AnalysisInteractive()
	analysis.cli()
