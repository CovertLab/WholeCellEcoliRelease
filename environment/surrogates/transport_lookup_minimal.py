from __future__ import absolute_import, division, print_function

import os
import csv
import time
import math
from scipy import constants

from agent.inner import CellSimulation
from reconstruction.spreadsheets import JsonReader

TUMBLE_JITTER = 2.0 # (radians)
DEFAULT_COLOR = [color/255 for color in [255, 69, 0]]

CSV_DIALECT = csv.excel_tab
TRANSPORT_REACTIONS_FILE = os.path.join("environment", "condition", "look_up_tables", "transport_reactions.tsv")
LIST_OF_LOOKUP_FILES = (
	os.path.join("environment", "condition", "look_up_tables", "avg_flux", "minimal.tsv"),
	os.path.join("environment", "condition", "look_up_tables", "avg_flux", "minimal_minus_oxygen.tsv"),
	os.path.join("environment", "condition", "look_up_tables", "avg_flux", "minimal_plus_amino_acids.tsv"),
)

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

class TransportMinimal(CellSimulation):
	''''''

	def __init__(self, state):
		self.initial_time = state.get('time', 0.0)
		self.local_time = state.get('time', 0.0)
		self.media_id = state.get('media_id', 'minimal')
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
		self.all_transport_reactions = {}
		with open(TRANSPORT_REACTIONS_FILE, 'rU') as csvfile:
			reader = JsonReader(csvfile, dialect=CSV_DIALECT)
			for row in reader:
				reaction_id = row["reaction id"]
				stoichiometry = row["stoichiometry"]
				reversible = row["is reversible"]
				catalyzed = row["catalyzed by"]
				self.all_transport_reactions[reaction_id] = {
					"stoichiometry": stoichiometry,
					"is reversible": reversible,
					"catalyzed by": catalyzed,
				}

		# make a dictionary with saved average fluxes for all transport reactions, in the three conditions
		# fluxes are in mmol/L
		self.flux_lookup = {}
		for file_name in LIST_OF_LOOKUP_FILES:
			attrName = file_name.split(os.path.sep)[-1].split(".")[0]
			self.flux_lookup[attrName] = {}
			with open(file_name, 'rU') as csvfile:
				reader = JsonReader(csvfile, dialect=CSV_DIALECT)
				for row in reader:
					reaction_id = row["reaction id"]
					flux = row["average flux mmol/L"]
					self.flux_lookup[attrName][reaction_id] = flux

		# exchange_ids declares which molecules' exchange will be controlled by transport
		aa_p_ids = [aa_id + "[p]" for aa_id in amino_acids]
		exchange_molecules = ["OXYGEN-MOLECULE[p]", "GLC[p]"]
		exchange_ids = exchange_molecules + aa_p_ids
		self.transport_reactions_ids = self.reactions_from_exchange(exchange_ids)

		# get the current flux lookup table, and set initial transport fluxes
		self.current_flux_lookup = self.flux_lookup[self.media_id]
		self.transport_fluxes = self.get_fluxes(self.current_flux_lookup, self.transport_reactions_ids)


	def update_state(self):
		# nAvogadro is in 1/mol --> convert to 1/mmol. volume is in fL --> convert to L
		self.molar_to_counts = (self.nAvogadro * 1e-3) * (self.volume * 1e-15)

		# get transport fluxes
		self.transport_fluxes = self.get_fluxes(self.current_flux_lookup, self.transport_reactions_ids)

		# convert to counts
		delta_counts = self.flux_to_counts(self.transport_fluxes)

		environment_deltas = {}
		for molecule in self.external_concentrations.keys():
			# TODO -- use external exchange map rather than (molecule + '[p]')
			molecule_p = molecule + '[p]'
			if molecule_p in delta_counts:
				environment_deltas[molecule] = delta_counts[molecule_p]

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

		# update lookup table
		self.current_flux_lookup = self.flux_lookup[self.media_id]

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


	## Flux-related functions
	def get_fluxes(self, flux_lookup, transport_reactions_ids):
		# TODO -- get reversible reactions, some fluxes are negative
		transport_fluxes = {
			transport_id: max(flux_lookup[transport_id], 0.0)
			for transport_id in transport_reactions_ids}
		# transport_fluxes = self.adjust_fluxes(transport_fluxes)
		return transport_fluxes

	def adjust_fluxes(self, transport_fluxes):
		'''adjust fluxes found by look up table'''

		added_flux = 0  # 1e-2 * (1 + math.sin(10 * self.local_time))
		adjusted_transport_fluxes = {
			transport_id: max(flux + added_flux, 0.0)
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
