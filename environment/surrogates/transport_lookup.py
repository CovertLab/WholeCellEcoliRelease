from __future__ import absolute_import, division, print_function

import time
import os
import csv
import math
from scipy import constants

from agent.inner import CellSimulation
from environment.condition.look_up_tables.look_up import LookUp
from reconstruction.spreadsheets import JsonReader
from itertools import ifilter

EXTERNAL_MOLECULES_FILE = os.path.join('environment', 'condition', 'environment_molecules.tsv')
REACTIONS_FILE = os.path.join("reconstruction", "ecoli", "flat", "reactions.tsv")
TRANSPORT_IDS_FILE = os.path.join("reconstruction", "ecoli", "flat", "transport_reactions.tsv")

TSV_DIALECT = csv.excel_tab
TUMBLE_JITTER = 2.0 # (radians)
DEFAULT_COLOR = [color/255 for color in [255, 69, 0]]

amino_acids = [
	'L-ALPHA-ALANINE',
	'ARG',
	'ASN',
	'L-ASPARTATE',
	'CYS',
	'GLT',
	'GLN',
	'GLY',
	'HIS',
	'ILE',
	'LEU',
	'LYS',
	'MET',
	'PHE',
	'PRO',
	'SER',
	'THR',
	'TRP',
	'TYR',
	'L-SELENOCYSTEINE',
	'VAL'
]

aa_p_ids = [aa_id + "[p]" for aa_id in amino_acids]
exchange_molecules = ["OXYGEN-MOLECULE[p]", "GLC[p]"]
exchange_ids = exchange_molecules + aa_p_ids

class TransportLookup(CellSimulation):
	''''''

	def __init__(self, state):
		self.initial_time = state.get('time', 0.0)
		self.local_time = state.get('time', 0.0)
		self.media_id = state.get('media_id', 'minimal')
		self.lookup_type = state.get('lookup', 'average')
		self.timestep = 1.0
		self.environment_change = {}
		self.volume = 1.0  # (fL)
		self.division_time = 100
		self.nAvogadro = constants.N_A

		# initial state
		self.external_concentrations = {}
		self.internal_concentrations = {}
		self.motile_force = [0.01, 0.01] # initial magnitude and relative orientation
		self.division = []

		# make dict of transport reactions
		# get all reactions
		all_reactions = {}
		with open(REACTIONS_FILE, 'rU') as tsvfile:
			reader = JsonReader(
				ifilter(lambda x: x.lstrip()[0] != "#", tsvfile), # Strip comments
				dialect = TSV_DIALECT)
			for row in reader:
				reaction_id = row["reaction id"]
				stoichiometry = row["stoichiometry"]
				reversible = row["is reversible"]
				catalyzed = row["catalyzed by"]
				all_reactions[reaction_id] = {
					"stoichiometry": stoichiometry,
					"is reversible": reversible,
					"catalyzed by": catalyzed,
				}

		# make dict of reactions in TRANSPORT_IDS_FILE
		self.all_transport_reactions = {}
		with open(TRANSPORT_IDS_FILE, 'rU') as tsvfile:
			reader = JsonReader(
				ifilter(lambda x: x.lstrip()[0] != "#", tsvfile), # Strip comments
				dialect = TSV_DIALECT)
			for row in reader:
				reaction_id = row["reaction id"]
				self.all_transport_reactions[reaction_id] = {
					"stoichiometry": all_reactions[reaction_id]["stoichiometry"],
					"is reversible": all_reactions[reaction_id]["is reversible"],
					"catalyzed by": all_reactions[reaction_id]["catalyzed by"],
				}

		# Make map of external molecule_ids with a location tag (as used in reaction stoichiometry) to molecule_ids in the environment
		self.molecule_to_external_map = {}
		self.external_to_molecule_map = {}
		with open(EXTERNAL_MOLECULES_FILE, 'rU') as tsvfile:
			reader = JsonReader(
				ifilter(lambda x: x.lstrip()[0] != "#", tsvfile), # Strip comments
				dialect = TSV_DIALECT)
			for row in reader:
				molecule_id = row['molecule id']
				location = row['exchange molecule location']
				self.molecule_to_external_map[molecule_id + location] = molecule_id
				self.external_to_molecule_map[molecule_id] = molecule_id + location

		# exchange_ids declares which molecules' exchange will be applied
		self.transport_reaction_ids = self.reactions_from_exchange(exchange_ids)

		# make look up object
		self.look_up = LookUp()

		# get the fluxes
		self.transport_fluxes = self.look_up.look_up(
			self.lookup_type,
			self.media_id,
			self.transport_reaction_ids)

		# adjust the fluxes
		# self.transport_fluxes = self.adjust_fluxes(self.transport_fluxes)

	def update_state(self):
		# nAvogadro is in 1/mol --> convert to 1/mmol. volume is in fL --> convert to L
		self.molar_to_counts = (self.nAvogadro * 1e-3) * (self.volume * 1e-15)

		# get transport fluxes
		self.transport_fluxes = self.look_up.look_up(
			self.lookup_type,
			self.media_id,
			self.transport_reaction_ids)
		# self.transport_fluxes = self.adjust_fluxes(self.transport_fluxes)

		# convert to counts
		delta_counts = self.flux_to_counts(self.transport_fluxes)

		# Get the deltas for environmental molecules
		environment_deltas = {}
		for molecule_id in delta_counts.keys():
			if molecule_id in self.molecule_to_external_map:
				external_molecule_id = self.molecule_to_external_map[molecule_id]
				environment_deltas[external_molecule_id] = delta_counts[molecule_id]

		# accumulate in environment_change
		self.accumulate_deltas(environment_deltas)

	def accumulate_deltas(self, environment_deltas):
		for molecule_id, count in environment_deltas.iteritems():
			self.environment_change[molecule_id] += count

	def check_division(self):
		# update division state based on time since initialization
		if self.local_time >= self.initial_time + self.division_time:
			self.division = [{'time': self.local_time}, {'time': self.local_time}]
		return self.division

	def time(self):
		return self.local_time

	def apply_outer_update(self, update):
		self.external_concentrations = update['concentrations']
		self.media_id = update['media_id']

		# reset environment change
		self.environment_change = {}
		for molecule in self.external_concentrations.iterkeys():
			self.environment_change[molecule] = 0

	def run_incremental(self, run_until):
		'''run until run_until'''
		while self.time() < run_until:
			self.local_time += self.timestep
			self.update_state()
			# self.check_division()

		time.sleep(1.0)  # pause for better coordination with Lens visualization. TODO: remove this

	def generate_inner_update(self):
		return {
			'volume': self.volume,
			'motile_force': self.motile_force,
			'environment_change': self.environment_change,
			'division': self.division,
			'color': DEFAULT_COLOR,
			'transport_fluxes': self.transport_fluxes,
			}

	def adjust_fluxes(self, transport_fluxes):
		'''adjust fluxes found by look up table'''

		time_constant = 10
		adjusted_transport_fluxes = {
			transport_id: max(flux(1 + math.sin(time_constant * self.local_time)), 0.0)
			for transport_id, flux in transport_fluxes.iteritems()}

		return adjusted_transport_fluxes

	def flux_to_counts(self, fluxes):
		rxn_counts = {
			reaction_id: int(self.molar_to_counts * flux)
			for reaction_id, flux in fluxes.iteritems()}
		delta_counts = {}
		for reaction_id, rxn_count in rxn_counts.iteritems():
			stoichiometry = self.all_transport_reactions[reaction_id]['stoichiometry']
			substrate_counts = {
				substrate_id: coeff * rxn_count
				for substrate_id, coeff in stoichiometry.iteritems()}
			# add to delta_counts
			for substrate, delta in substrate_counts.iteritems():
				if substrate in delta_counts:
					delta_counts[substrate] += delta
				else:
					delta_counts[substrate] = delta
		return delta_counts

	def reactions_from_exchange(self, include_exchanges):
		include_reactions = []
		for reaction_id, specs in self.all_transport_reactions.iteritems():
			reaction_molecules = specs['stoichiometry'].keys()
			for exchange in include_exchanges:
				if exchange in reaction_molecules:
					include_reactions.append(reaction_id)
		return include_reactions
