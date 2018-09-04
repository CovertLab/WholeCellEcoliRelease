from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

import os

import numpy as np
from scipy import constants

in_sherlock = 'SHERLOCK' in os.environ

import matplotlib

# Turn off interactive plotting when running on sherlock
if not in_sherlock:
	matplotlib.use('TKAgg')

import matplotlib.pyplot as plt

if not in_sherlock:
	plt.ion()
	fig = plt.figure()

# Constants
N_AVOGADRO = constants.N_A # TODO (ERAN) get this from sim_data.constants.nAvogadro
PI = np.pi

# Lattice parameters
N_DIMS = 2
PATCHES_PER_EDGE = 10
TOTAL_VOLUME = 1E-11  # (L)
EDGE_LENGTH = 10.  # (micrometers). for reference: e.coli length is on average 2 micrometers.

# Physical constants
DIFFUSION = 0.001  # diffusion constant. (micrometers^2/s) # TODO (Eran) calculate correct diffusion rate

# Derived environmental constants
PATCH_VOLUME = TOTAL_VOLUME / (PATCHES_PER_EDGE*PATCHES_PER_EDGE)
DX = EDGE_LENGTH / PATCHES_PER_EDGE  # intervals in x- directions (assume y- direction equivalent)
DX2 = DX*DX
# DT = DX2 * DX2 / (2 * DIFFUSION * (DX2 + DX2)) # upper limit on the time scale (go with at least 50% of this)

# Cell constants
CELL_RADIUS = 0.5 # (micrometers)
ORIENTATION_JITTER = PI/40  # (radians/s)
LOCATION_JITTER = 0.01 # (micrometers/s)

class EnvironmentSpatialLattice(object):
	def __init__(self, concentrations):
		self._time = 0
		self._timestep = 1.0
		self.run_for = 5

		self.agent_id = -1

		self.simulations = {}
		self.locations = {}
		self.volumes = {}

		self._molecule_ids = concentrations.keys()
		self.concentrations = concentrations.values()

		# Create lattice and fill each site with concentrations dictionary
		# Molecule identities are defined along the major axis, with spatial dimensions along the other two axes.
		self.lattice = np.empty([len(self._molecule_ids)] + [PATCHES_PER_EDGE for dim in xrange(N_DIMS)], dtype=np.float64)
		for idx, molecule in enumerate(self._molecule_ids):
			self.lattice[idx].fill(self.concentrations[idx])

		if os.path.exists("out/manual/environment.txt"):
			os.remove("out/manual/environment.txt")
		if os.path.exists("out/manual/locations.txt"):
			os.remove("out/manual/locations.txt")

		glucose_lattice = self.lattice[self._molecule_ids.index('GLC[p]')]
		plt.imshow(glucose_lattice, vmin=0, vmax=25, cmap='YlGn')
		plt.colorbar()
		plt.axis('off')
		plt.pause(0.0001)


	def evolve(self):
		''' Evolve environment '''

		# TODO (Eran) updating location with the environment causes cells to jump at run_for ts, and then cells deposit
		# all of their deltas at the jumps, rather than along their path laid out during evolve()
		self.update_locations()

		self.run_diffusion()


	def update_locations(self):
		''' Update location for all agent_ids '''
		for agent_id, location in self.locations.iteritems():
			# Move the cell around randomly
			self.locations[agent_id][0:2] = (location[0:2] +
				np.random.normal(scale=np.sqrt(LOCATION_JITTER * self._timestep), size=N_DIMS)
				) % EDGE_LENGTH

			# Orientation jitter
			self.locations[agent_id][2] = (location[2] +
				np.random.normal(scale=ORIENTATION_JITTER * self._timestep)
				)

	def run_diffusion(self):
		change_lattice = np.zeros(self.lattice.shape)
		for idx in xrange(len(self.lattice)):
			molecule = self.lattice[idx]

			# run diffusion if molecule field is not uniform
			if (len(set(molecule.flatten())) != 1):
				change_lattice[idx] = self.diffusion_timestep(molecule)

		self.lattice += change_lattice


	def diffusion_timestep(self, lattice):
		''' calculate concentration changes cause by diffusion. Assumes periodic lattice, with wrapping'''

		# TODO (Eran) write this as matrix operation rather than np.roll.
		N = np.roll(lattice, 1, axis=0)
		S = np.roll(lattice, -1, axis=0)
		W = np.roll(lattice, 1, axis=1)
		E = np.roll(lattice, -1, axis=1)

		change_lattice = DIFFUSION * self._timestep * ((N + S + W + E - 4 * lattice) / DX2)

		return change_lattice


	def run_incremental(self, run_until):
		''' Simulate until run_until '''
		self.output_environment()
		self.output_locations()

		while self._time < run_until:
			self._time += self._timestep
			self.evolve()


	def output_environment(self):
		'''plot environment lattice'''
		glucose_lattice = self.lattice[self._molecule_ids.index('GLC[p]')]

		plt.clf()
		plt.imshow(glucose_lattice, cmap='YlGn')
		plt.colorbar()
		plt.axis('off')


	def output_locations(self):
		'''plot cell locations and orientations'''
		for agent_id, location in self.locations.iteritems():
			y = location[0] * PATCHES_PER_EDGE / EDGE_LENGTH - 0.5
			x = location[1] * PATCHES_PER_EDGE / EDGE_LENGTH - 0.5
			theta = location[2]
			volume = self.volumes[agent_id]

			# get length, scaled to lattice resolution
			length = self.volume_to_length(volume)

			dx = length * PATCHES_PER_EDGE / EDGE_LENGTH * np.sin(theta)
			dy = length * PATCHES_PER_EDGE / EDGE_LENGTH * np.cos(theta)

			plt.plot([x-dx/2, x+dx/2], [y-dy/2, y+dy/2],
				color='salmon', linewidth=CELL_RADIUS/EDGE_LENGTH*600, solid_capstyle='round')
			plt.plot([x-dx*9.5/20, x+dx*9.5/20], [y-dy*9.5/20, y+dy*9.5/20],
				color='slateblue', linewidth=CELL_RADIUS/EDGE_LENGTH*450, solid_capstyle='round')

		if not in_sherlock:
			plt.pause(0.0001)


	def volume_to_length(self, volume):
		'''
		get cell length from volume, using the following equation for capsule volume, with V=volume, r=radius,
		a=length of cylinder without rounded caps, l=total length:

		V = (4/3)*PI*r^3 + PI*r^2*a
		l = a + 2*r
		'''

		cylinder_length = (volume - (4/3) * PI * CELL_RADIUS**3) / (PI * CELL_RADIUS**2)
		total_length = cylinder_length + 2 * CELL_RADIUS

		return total_length


	def counts_to_concentration(self, counts):
		''' Convert an array of counts to concentrations '''
		concentrations = [count / (PATCH_VOLUME * N_AVOGADRO) for count in counts]
		return concentrations


	def update_from_simulations(self, all_changes):
		'''
		Use change counts from all the inner simulations, convert them to concentrations,
		and add to the environmental concentrations of each molecule at each simulation's location
		'''
		for agent_id, update in all_changes.iteritems():
			self.volumes[agent_id] = update['volume']
			change_counts = update['environment_change']

			location = self.locations[agent_id][0:2] * PATCHES_PER_EDGE / EDGE_LENGTH
			patch_site = tuple(np.floor(location).astype(int))

			change_concentrations = self.counts_to_concentration(change_counts.values())
			for molecule, change_conc in zip(change_counts.keys(), change_concentrations):
				self.lattice[self._molecule_ids.index(molecule), patch_site[0], patch_site[1]] += change_conc


	def get_molecule_ids(self):
		''' Return the ids of all molecule species in the environment '''
		return self._molecule_ids


	def get_concentrations(self):
		'''returns a dict with {molecule_id: conc} for each sim give its current location'''
		concentrations = {}
		for agent_id in self.simulations.keys():
			# get concentration from cell's given bin
			location = self.locations[agent_id][0:2] * PATCHES_PER_EDGE / EDGE_LENGTH
			patch_site = tuple(np.floor(location).astype(int))

			concentrations[agent_id] = dict(zip(self._molecule_ids, self.lattice[:,patch_site[0],patch_site[1]]))
		return concentrations


	def time(self):
		return self._time


	def add_simulation(self, agent_id):
		state = {}

		# Place cell at a random initial location
		location = np.random.uniform(0,EDGE_LENGTH,N_DIMS)
		orientation = np.random.uniform(0, 2*PI)

		self.simulations[agent_id] = state
		self.locations[agent_id] = np.hstack((location, orientation))
		self.volumes[agent_id] = 1.


	def remove_simulation(self, agent_id):
		self.simulations.pop(agent_id, {})
		self.locations.pop(agent_id, {})
		self.volumes.pop(agent_id, {})


	def run_simulations_until(self):
		until = {}
		run_until = self.time() + self.run_for
		for agent_id in self.simulations.keys():
			until[agent_id] = run_until

		# Pass the environment a run_until
		until[self.agent_id] = run_until

		return until
