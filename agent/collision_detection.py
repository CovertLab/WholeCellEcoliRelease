from __future__ import absolute_import, division, print_function

import numpy as np

import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
plt.ion()
fig = plt.figure()
from agent.grid import Grid, Rectangle

def collision_detection(agents, lattice_size, dx):
	grid = Grid([lattice_size, lattice_size], dx)

	for agent_id, agent in agents.iteritems():
		box = agent['indices']
		grid.impress(box)

	total_overlap = grid.overlap()

	# get forces
	forces = {}
	for agent_id, agent in agents.iteritems():
		location = agent['location']
		box = agent['indices']
		forces[agent_id] = grid.get_forces(location, box)

	return grid, total_overlap, forces


def accept(delta, temp):
	prob_accept = np.exp(-delta/temp)

	return prob_accept


def volume_exclusion(agents):
	introspect = {}

	grid, overlap, forces = collision_detection(agents, edge_length, resolution)
	cycles = 0
	while overlap > 0:
		agents_new = agents.copy()
		introspect['searches'] = {}

		# update one agent at a time
		for agent_id, specs in agents_new.iteritems():
			agent = agents_new[agent_id]
			force = forces[agent_id]

			location = agent['location'] + force + np.random.normal(scale=np.sqrt(TRANSLATIONAL_JITTER), size=2)
			orientation = agent['orientation'] + np.random.normal(scale=ROTATIONAL_JITTER) % (2 * np.pi)

			# check if agent is in the environmental bounds
			box = Rectangle([agent['radius'], agent['length']], location, orientation)
			indices = box.render(resolution)
			in_bounds = grid.check_in_bounds(indices)

			if in_bounds:
				agent['location'] = location
				agent['orientation'] = orientation
				agent['indices'] = indices

		grid_new, overlap_new, forces_new = collision_detection(agents_new, edge_length, resolution)

		if overlap_new <= overlap:
			agents = agents_new
			grid = grid_new
			overlap = overlap_new
			forces = forces_new

			plt.imshow(grid.grid)
			plt.pause(0.0001)

		elif np.random.rand() < accept(overlap_new - overlap, TEMPERATURE):
			agents = agents_new
			grid = grid_new
			overlap = overlap_new
			forces = forces_new

			plt.imshow(grid.grid)
			plt.pause(0.0001)


ROTATIONAL_JITTER = 0.01 # (radians/s)
TRANSLATIONAL_JITTER = 0.001 # (micrometers/s)
TEMPERATURE = 20 # for acceptance function

edge_length = 10
resolution = 0.1
agents = {}
animals = [
	'aardvark',
	'basilisk',
	'capybara',
	'dingo',
	'elephant',
	'fox',
	'groundhog',
	'hyena',
	'ibis',
	'jackal',
	'koala',
	'llama',
	'marmuset',
	'narwhal',
	'ocelot']
	
for animal in animals:
	agents[animal] = {
		'location': np.random.random(2) * 7 + 1.5,
		'orientation': np.random.random(1)[0] * np.pi * 2,
		'length': 2,
		'radius': 0.5}

#give agents indices
for agent_id, agent in agents.iteritems():
	box = Rectangle([agent['radius'], agent['length']], agent['location'], agent['orientation'])
	indices = box.render(resolution)
	agent['indices'] = indices

volume_exclusion(agents)


import ipdb; ipdb.set_trace()
