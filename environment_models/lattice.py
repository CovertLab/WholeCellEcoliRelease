"""
Lattice

A two-dimensional lattice environmental model

## Physics
# Diffusion constant of glucose in 0.5 and 1.5 percent agarose gel = ~6 * 10^-10 m^2/s (Weng et al. 2005. Transport of glucose and poly(ethylene glycol)s in agarose gels).
# Conversion to micrometers: 6 * 10^-10 m^2/s = 600 micrometers^2/s.

## Cell biophysics
# rotational diffusion in liquid medium with viscosity = 1 mPa.s: Dr = 3.5+/-0.3 rad^2/s (Saragosti, et al. 2012. Modeling E. coli tumbles by rotational diffusion.)
# translational diffusion in liquid medium with viscosity = 1 mPa.s: Dt=100 micrometers^2/s (Saragosti, et al. 2012. Modeling E. coli tumbles by rotational diffusion.)


@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

import os

import numpy as np
from scipy import constants

animating = 'ENVIRONMENT_ANIMATION' in os.environ

import matplotlib

# Turn off interactive plotting when running on sherlock
if animating:
	matplotlib.use('TKAgg')

import matplotlib.pyplot as plt

from agent.outer import EnvironmentSimulation


if animating:
	plt.ion()
	fig = plt.figure()

# Constants
N_AVOGADRO = constants.N_A
PI = np.pi

# Lattice parameters
N_DIMS = 2
PATCHES_PER_EDGE = 10 # TODO (Eran) this should scale to accomodate diffusion

EDGE_LENGTH = 10.0  # (micrometers)
DEPTH = 3000.0 # (micrometers). An average Petri dish has a depth of 3-4 mm
TOTAL_VOLUME = (DEPTH * EDGE_LENGTH**2) * (10**-15) # (L)

# Physical constants
DIFFUSION = 0.1  # (micrometers^2/s)

# Derived environmental constants
PATCH_VOLUME = TOTAL_VOLUME / (PATCHES_PER_EDGE**2)
DX = EDGE_LENGTH / PATCHES_PER_EDGE  # intervals in x- directions (assume y- direction equivalent)
DX2 = DX*DX
# DT = DX2 * DX2 / (2 * DIFFUSION * (DX2 + DX2)) # upper limit on the time scale (go with at least 50% of this)

# Cell constants
CELL_RADIUS = 0.5 # (micrometers)
ROTATIONAL_JITTER = 0.05 # (radians/s)
TRANSLATIONAL_JITTER = 0.001 # (micrometers/s)

class EnvironmentSpatialLattice(EnvironmentSimulation):
	def __init__(self, concentrations):
		self._time = 0
		self._timestep = 1.0
		self._run_for = 5

		self.simulations = {}  # map of agent_id to simulation state
		self.locations = {}    # map of agent_id to location and orientation

		self._molecule_ids = concentrations.keys()
		self.concentrations = concentrations.values()
		self.molecule_index = {molecule: index for index, molecule in enumerate(self._molecule_ids)}

		# Create lattice and fill each site with concentrations dictionary
		# Molecule identities are defined along the major axis, with spatial dimensions along the other two axes.
		self.lattice = np.empty([len(self._molecule_ids)] + [PATCHES_PER_EDGE for dim in xrange(N_DIMS)], dtype=np.float64)
		for idx, molecule in enumerate(self._molecule_ids):
			self.lattice[idx].fill(self.concentrations[idx])

		if os.path.exists("out/manual/environment.txt"):
			os.remove("out/manual/environment.txt")
		if os.path.exists("out/manual/locations.txt"):
			os.remove("out/manual/locations.txt")

		if animating:
			glucose_lattice = self.lattice[self.molecule_index['GLC[p]']]
			plt.imshow(glucose_lattice, vmin=0, vmax=25, cmap='YlGn')
			plt.colorbar()
			plt.axis('off')
			plt.pause(0.0001)

	def evolve(self):
		''' Evolve environment '''

		self.update_locations()

		self.run_diffusion()


	def update_locations(self):
		''' Update location for all agent_ids '''
		for agent_id, location in self.locations.iteritems():

			# Translational diffusion
			self.locations[agent_id][0:2] += np.random.normal(scale=np.sqrt(TRANSLATIONAL_JITTER * self._timestep), size=N_DIMS)

			# Bounce cells off of lattice edges
			self.locations[agent_id][0:2][self.locations[agent_id][0:2] >= EDGE_LENGTH] -= 2 * self.locations[agent_id][0:2][self.locations[agent_id][0:2]>= EDGE_LENGTH] % EDGE_LENGTH

			# Rotational diffusion
			self.locations[agent_id][2] = (location[2] + np.random.normal(scale=ROTATIONAL_JITTER * self._timestep)) % (2 * PI)


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
		if animating:
			self.output_environment()
			self.output_locations()

		while self._time < run_until:
			self._time += self._timestep
			self.evolve()


	def output_environment(self):
		'''plot environment lattice'''
		glucose_lattice = self.lattice[self.molecule_index['GLC[p]']]

		plt.clf()
		# plt.imshow(np.pad(glucose_lattice, ((1,1),(1,1)), 'wrap'), cmap='YlGn')
		plt.imshow(glucose_lattice, cmap='YlGn')
		plt.title('time: ' + str(self._time) + ' (s)')
		plt.colorbar()
		plt.axis('off')
		# plt.ylim((-DX, EDGE_LENGTH+DX))
		# plt.xlim((-DX, EDGE_LENGTH+DX))


	def output_locations(self):
		'''plot cell locations and orientations'''
		for agent_id, location in self.locations.iteritems():
			y = location[0] * PATCHES_PER_EDGE / EDGE_LENGTH
			x = location[1] * PATCHES_PER_EDGE / EDGE_LENGTH
			theta = location[2]
			volume = self.simulations[agent_id]['state']['volume']

			# get length, scaled to lattice resolution
			length = self.volume_to_length(volume)

			dx = length * PATCHES_PER_EDGE / EDGE_LENGTH * np.sin(theta)
			dy = length * PATCHES_PER_EDGE / EDGE_LENGTH * np.cos(theta)

			plt.plot([x-dx/2, x+dx/2], [y-dy/2, y+dy/2],
				color='slateblue', linewidth=CELL_RADIUS/EDGE_LENGTH*600, solid_capstyle='round')


		if animating:
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


	def count_to_concentration(self, count):
		''' Convert count to concentrations '''
		return count / (PATCH_VOLUME * N_AVOGADRO)


	def apply_inner_update(self, update, now):
		'''
		Use change counts from all the inner simulations, convert them to concentrations,
		and add to the environmental concentrations of each molecule at each simulation's location
		'''
		self.simulations.update(update)

		for agent_id, simulation in self.simulations.iteritems():
			# only apply changes if we have reached this simulation's time point.
			if simulation['time'] <= now:
				print('=================== simulation update: {}'.format(simulation))
				state = simulation['state']
				location = self.locations[agent_id][0:2] * PATCHES_PER_EDGE / EDGE_LENGTH
				patch_site = tuple(np.floor(location).astype(int))

				for molecule, count in state['environment_change'].iteritems():
					concentration = self.count_to_concentration(count)
					index = self.molecule_index[molecule]
					self.lattice[index, patch_site[0], patch_site[1]] += concentration


	def get_molecule_ids(self):
		''' Return the ids of all molecule species in the environment '''
		return self._molecule_ids


	def generate_outer_update(self, now):
		'''returns a dict with {molecule_id: conc} for each sim give its current location'''
		update = {}
		for agent_id, simulation in self.simulations.iteritems():
			# only provide concentrations if we have reached this simulation's time point.
			if simulation['time'] <= now:
				# get concentration from cell's given bin
				location = self.locations[agent_id][0:2] * PATCHES_PER_EDGE / EDGE_LENGTH
				patch_site = tuple(np.floor(location).astype(int))
				update[agent_id] = {}
				update[agent_id]['concentrations'] = dict(zip(
					self._molecule_ids,
					self.lattice[:,patch_site[0],patch_site[1]]))

		return update


	def time(self):
		return self._time

	def add_simulation(self, agent_id, simulation):
		if agent_id not in self.simulations:
			self.simulations[agent_id] = {}
		self.simulations[agent_id].update(simulation)

		if agent_id not in self.locations:
			# Place cell at either the provided or a random initial location
			location = simulation['agent_config'].get(
				'location', np.random.uniform(0,EDGE_LENGTH,N_DIMS))
			orientation = simulation['agent_config'].get(
				'orientation', np.random.uniform(0, 2*PI))

			self.locations[agent_id] = np.hstack((location, orientation))

	def simulation_parameters(self, agent_id):
		latest = max([
			simulation['time']
			for agent_id, simulation
			in self.simulations.iteritems()])
		time = max(self._time, latest)
		return {'time': time}

	def simulation_state(self, agent_id):
		return dict(
			self.simulations[agent_id],
			location=self.locations[agent_id])

	def rotation_matrix(self, orientation):
		sin = np.sin(orientation)
		cos = np.cos(orientation)
		return np.matrix([
			[cos, -sin],
			[sin, cos]])

	def daughter_location(self, location, orientation, length, index):
		offset = np.array([length * 0.75, 0])
		rotation = self.rotation_matrix(-orientation + (index * np.pi))
		translation = (offset * rotation).A1
		return (location + translation)

	def apply_parent_state(self, agent_id, parent):
		parent_location = parent['location']
		index = parent['index']
		orientation = parent_location[2]
		volume = self.simulations[agent_id]['state']['volume'] * 0.5
		length = self.volume_to_length(volume)
		location = self.daughter_location(parent_location[0:2], orientation, length, index)

		print("================= parent: {} - daughter: {}".format(parent_location, location))

		self.locations[agent_id] = np.hstack((location, orientation))

	def run_for(self):
		return self._run_for

	def remove_simulation(self, agent_id):
		self.simulations.pop(agent_id, {})
		self.locations.pop(agent_id, {})
