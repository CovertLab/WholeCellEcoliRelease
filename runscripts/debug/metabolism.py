#! /usr/bin/env python
"""
Tools and analysis to debug metabolism problems.
"""

from __future__ import absolute_import, division, print_function

import argparse
import os
from typing import List

import numpy as np
from six.moves import cPickle, range, zip

from models.ecoli.processes.metabolism import (CONC_UNITS, CONVERSION_UNITS,
	FluxBalanceAnalysisModel, GDCW_BASIS)
from wholecell.io.tablereader import TableReader
from wholecell.utils import constants, scriptBase, units


class MetabolismDebug(scriptBase.ScriptBase):
	def define_parameters(self, parser):
		super(MetabolismDebug, self).define_parameters(parser)

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

		# Debug options
		parser.add_argument('--validation', type=int, default=1,
			help='Number of time steps to run for validation. If < 0, will run all.')
		parser.add_argument('-t', '--timestep', type=int, default=0,
			help='Specify a timestep to run analysis for.  Works with certain options.')
		parser.add_argument('--produce-met', nargs='*',
			help='Check which enzymes need to be positive to produce a list of metabolites.')

	def update_args(self, args):
		super(MetabolismDebug, self).update_args(args)

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

	def load_data(self, sim_data_file, sim_out_dir):
		# type: (str, str) -> None
		"""
		Loads sim_data and simulation output data from files and saves it as
		instance variables.

		Args:
			sim_data_file: path to the sim_data file for the simulation
			sim_out_dir: path to the simOut dir for the simulation
		"""

		with open(sim_data_file, 'rb') as f:
			self.sim_data = cPickle.load(f)
		self.exchange_molecules = np.array(self.sim_data.external_state.all_external_exchange_molecules)

		# Listeners used
		fba_reader = TableReader(os.path.join(sim_out_dir, 'FBAResults'))
		kinetics_reader = TableReader(os.path.join(sim_out_dir, 'EnzymeKinetics'))
		main_reader = TableReader(os.path.join(sim_out_dir, 'Main'))

		# Load data
		self.catalyst_ids = fba_reader.readAttribute('catalyst_ids')
		self.update_molecules = fba_reader.readAttribute('conc_update_molecules')
		self.media_ids = fba_reader.readColumn('media_id')[1:]
		self.conc_updates = fba_reader.readColumn('conc_updates')[1:, :]
		self.catalyst_counts = fba_reader.readColumn('catalyst_counts')[1:, :]
		self.translation_gtp = fba_reader.readColumn('translation_gtp')[1:]
		self.coefficients = CONVERSION_UNITS * fba_reader.readColumn('coefficient')[1:]
		self.unconstrained_molecules = fba_reader.readColumn('unconstrained_molecules')[1:, :]
		self.constrained_molecules = fba_reader.readColumn('constrained_molecules')[1:, :]
		self.uptake_constraints = fba_reader.readColumn('uptake_constraints')[1:, :]
		self.objective = fba_reader.readColumn('objectiveValue')[1:]

		self.metabolite_ids = kinetics_reader.readAttribute('metaboliteNames')
		self.metabolite_counts = kinetics_reader.readColumn('metaboliteCountsInit')[1:, :]
		self.counts_to_molar = CONC_UNITS * kinetics_reader.readColumn('countsToMolar')[1:]

		self.time_step_sizes = units.s * main_reader.readColumn('timeStepSec')[1:]
		self.n_time_steps = len(self.time_step_sizes)

	def new_model(self):
		"""Create new model."""
		return FluxBalanceAnalysisModel(self.sim_data)

	def solve_timestep(self, model, timestep):
		# type: (FluxBalanceAnalysisModel, int) -> None
		"""
		Solves the FBA problem for a given timestep.

		Args:
			model: FBA model to solve
			timestep: simulation timestep to select data from
		"""

		# Calculations
		conc_updates = dict(zip(self.update_molecules, self.conc_updates[timestep, :]))
		catalyst_dict = dict(zip(self.catalyst_ids, self.catalyst_counts[timestep, :]))
		metabolite_dict = dict(zip(self.metabolite_ids, self.metabolite_counts[timestep, :]))
		unconstrained = set(self.exchange_molecules[self.unconstrained_molecules[timestep, :]])
		constrained = {
			mol: uptake * GDCW_BASIS
			for mol, uptake, present in
			zip(self.exchange_molecules, self.uptake_constraints[timestep, :], self.constrained_molecules[timestep, :])
			if present
			}
		kinetic_enzyme_counts = np.array([catalyst_dict[e]
			for e in model.kinetic_constraint_enzymes])
		kinetic_substrate_counts = np.array([metabolite_dict[s]
			for s in model.kinetic_constraint_substrates])
		counts_to_molar = self.counts_to_molar[timestep]
		coefficient = self.coefficients[timestep]

		# Set molecule availability (internal and external)
		model.set_molecule_levels(
			self.metabolite_counts[timestep, :], counts_to_molar,
			coefficient, self.media_ids[timestep],
			unconstrained, constrained, conc_updates,
			)

		# Set reaction limits for maintenance and catalysts present
		model.set_reaction_bounds(
			self.catalyst_counts[timestep, :], counts_to_molar,
			coefficient, self.translation_gtp[timestep],
			)

		# Constrain reactions based on targets
		model.set_reaction_targets(
			kinetic_enzyme_counts, kinetic_substrate_counts,
			counts_to_molar, self.time_step_sizes[timestep],
			)

		# Solve multiple times (if needed) for robustness
		n_retries = 3
		model.fba.solve(n_retries)

	def validation(self, n_steps):
		# type: (int) -> None
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
		model = self.new_model()
		for timestep in range(n_steps):
			self.solve_timestep(model, timestep)
			if model.fba.getObjectiveValue() != self.objective[timestep]:
				raise ValueError('Objective value does not match for time step {}'.format(timestep))
		print('All {} timesteps match the results from the whole-cell model.'.format(n_steps))

	def produce_metabolites(self, metabolites, timestep):
		# type: (List[str], int) -> None
		"""
		Check if metabolites can be produced by providing a missing enzyme.

		Args:
			metabolites: metabolite IDs to check (with a location tag)
			timestep: the simulation timestep to run

		TODO:
			- check for all metabolites that are initially not produced
			- check for enzyme combinations (eg need 2 instead of just 1 for flux)
		"""

		# Setup model to use
		model = self.new_model()
		metabolites_indices = {
			metabolite: idx
			for idx, metabolite in enumerate(model.fba.getOutputMoleculeIDs())
			}

		# Check each metabolite to see if it can be produced if enzymes exist
		old_catalyst_counts = self.catalyst_counts[timestep, :].copy()
		for met in metabolites:
			# Skip metabolites without concentrations
			if met in metabolites_indices:
				met_idx = metabolites_indices[met]
			else:
				print('FBA does not solve for concentration of {}. Skipping.'.format(met))
				continue

			# Skip metabolites that are already produced
			self.solve_timestep(model, timestep)
			if model.fba.getOutputMoleculeLevelsChange()[met_idx] > 0:
				print('Positive flux already exists for {}. Skipping.'.format(met))
				continue

			# Check if positive counts (1) for each enzyme can produce the metabolite
			print('Checking for production of {}:'.format(met))
			for i, catalyst in enumerate(self.catalyst_ids):
				# Skip enzymes that already exist
				if old_catalyst_counts[i] != 0:
					continue

				# Set counts from 0 to 1
				self.catalyst_counts[timestep, i] = 1

				# Solve FBA problem and check if metabolite counts increase
				self.solve_timestep(model, timestep)
				if model.fba.getOutputMoleculeLevelsChange()[met_idx] > 0:
					print('\t{} causes positive flux'.format(catalyst))

				# Reset change
				self.catalyst_counts[timestep, :] = old_catalyst_counts

	def run(self, args):
		# type: (argparse.Namespace) -> None

		self.load_data(args.sim_data_file, args.sim_out_dir)
		self.validation(args.validation)
		if args.produce_met:
			self.produce_metabolites(args.produce_met, args.timestep)


if __name__ == '__main__':
	analysis = MetabolismDebug()
	analysis.cli()
