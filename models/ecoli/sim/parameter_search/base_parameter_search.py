"""
Base implementation for parameter search methods.  These methods should define
the parameters to search over and the objective to evaluate progress.  Specific
implementations should subclass from BaseParameterSearch and implement the
required functions.
"""

import abc
import os
import pickle
from typing import Any, Dict, Tuple, Union

import numpy as np

from wholecell.io.tablereader import TableReader
from wholecell.sim.simulation import DEFAULT_SIMULATION_KWARGS


DEFAULT_CLI_KWARGS = {
	'init_sims': 1,
	'generations': 1,
	}


class Parameter(abc.ABC):
	@abc.abstractmethod
	def get_param(self, data, default=None):
		pass

	@abc.abstractmethod
	def set_param(self, data, value):
		pass

	def __str__(self):
		return self._name


class RawParameter(Parameter):
	def __init__(self, attr: str, id_columns: Dict[str, Any], columns: Any, name: str):
		"""
		Create an object that can access and modify a desired property in raw_data.

		Args:
			attr: attribute of raw_data containing parameter to modify
			id_columns: column names (key) and values (value) to identify the row to modify
				eg. {'Gene': 'TrpR'} will modify the first row that has TrpR listed in column Gene
			columns: columns within the row to get or modify with several options for indexing
				str: will get/modify the value of the column directly
				List[str]: will get the average of the columns specified and set the same updated
					value for each one
				List[List[Any]]: will continue indexing into the column with each value in the inner
					list and get the average of the columns specified and set the same updated value
					for each one if multiple inner lists are provided (eg. [['KM', 0]] will change
					the first element in a list of KM values under column KM)
			name: the text identifier of this parameter to display for updates
		"""

		self._attr = attr
		self._id = id_columns
		self._columns = columns
		self._name = name
		self._cached_index = None

	def get_row(self, raw_data):
		def find_row(attr):
			for i, row in enumerate(attr):
				for col, val in self._id.items():
					if row[col] != val:
						break
				else:
					self._cached_index = i
					return row

			raise RuntimeError('Could not find a row matching the given ID columns.')

		attr = getattr(raw_data, self._attr)

		if self._cached_index is not None:
			row = attr[self._cached_index]
			for col, val in self._id.items():
				if row[col] != val:
					row = find_row(attr)
					break
		else:
			row = find_row(attr)

		return row

	def get_param(self, raw_data, default=None):
		row = self.get_row(raw_data)
		try:
			if isinstance(self._columns, list):
				values = []
				for column in self._columns:
					if isinstance(column, list):
						obj = row
						for subcolumn in column:
							obj = obj[subcolumn]
						values.append(obj)
					else:
						values.append(row[column])
				value = np.mean(values)
			else:
				value = row[self._columns]
		except (KeyError, IndexError):
			return default

		return value

	def set_param(self, raw_data, value):
		row = self.get_row(raw_data)
		if isinstance(self._columns, list):
			for column in self._columns:
				if isinstance(column, list):
					obj = row
					for subcolumn in column[:-1]:
						obj = obj[subcolumn]
					obj[column[-1]] = value
				else:
					row[column] = value
		else:
			row[self._columns] = value


class SimParameter(Parameter):
	def __init__(self, attr: str):
		"""
		Create an object that can access and modify a desired property in sim_data.

		Args:
			attr: attribute of sim_data to modify, dot separated for nested attributes
		"""

		self.attrs = attr.split('.')
		self._name = attr

	def get_param(self, obj, default=None):
		for a in self.attrs:
			if hasattr(obj, a):
				obj = getattr(obj, a)
			else:
				return default

		return obj

	def set_param(self, obj, value):
		for a in self.attrs[:-1]:
			obj = getattr(obj, a)
		setattr(obj, self.attrs[-1], value)


class BaseParameterSearch(abc.ABC):
	parca_args = {}  # type: Dict
	# TODO: handle raw and sim params the same - create a class for SimParameter and combine attributes below
	_raw_params = ()  # type: Union[Tuple, Tuple[RawParameter]]
	_sim_params = ()  # type: Union[Tuple, Tuple[SimParameter]]
	_init_raw_params = {}  # type: Dict
	_init_sim_params = {}  # type: Dict
	sims_to_run = ()  # type: Union[Tuple, Tuple[Dict]]

	def __init__(self):
		self.variant_name = self.__class__.__name__
		self.n_parameters = len(self._raw_params) + len(self._sim_params)
		self.raw_params = {p: None for p in self._raw_params}
		self.sim_params = {p: None for p in self._sim_params}
		self.initialized = False

	@abc.abstractmethod
	def get_objective(self, sim_out_dirs, sim_data_files):
		pass

	def initialize(self, raw_data_file, sim_data_file, iteration):
		# If no raw params, raw_data will not be saved with each iteration so
		# this would cause issues loading from a specific iteration
		if self.raw_params:
			with open(raw_data_file, 'rb') as f:
				raw_data = pickle.load(f)
			for param in self.raw_params:
				value = param.get_param(raw_data)
				if iteration == 0:
					self.raw_params[param] = self._init_raw_params.get(str(param), value)
				else:
					self.raw_params[param] = value

		with open(sim_data_file, 'rb') as f:
			sim_data = pickle.load(f)
		for param in self.sim_params:
			value = param.get_param(sim_data)
			if iteration == 0:
				self.sim_params[param] = self._init_sim_params.get(str(param), value)
			else:
				self.sim_params[param] = value

		self.initialized = True

	def get_sim_params(self, sim_dir, variants):
		all_params = []
		for variant in variants:
			for index, sim_params in enumerate(self.sims_to_run):
				params = DEFAULT_SIMULATION_KWARGS.copy()
				params.update(DEFAULT_CLI_KWARGS)
				params.update(sim_params)

				params['variant directory'] = os.path.join(sim_dir, f'{self.variant_name}_{variant:06n}')
				params['index'] = index

				all_params.append(params)

		return all_params

	def print_update(self):
		def print_params(params, label):
			if params:
				print(f'{label} parameters:')
				for p, val in params.items():
					print(f'\t{p}: {val}')

		print_params(self.raw_params, 'Raw data')
		print_params(self.sim_params, 'Sim data')

	def read_column(self, out_dir, table, column):
		return TableReader(os.path.join(out_dir, table)).readColumn(column)

	def read_attribute(self, out_dir, table, attr):
		return TableReader(os.path.join(out_dir, table)).readAttribute(attr)
