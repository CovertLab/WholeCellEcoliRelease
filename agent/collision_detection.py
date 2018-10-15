from __future__ import absolute_import, division, print_function

import numpy as np

import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
plt.ion()
fig = plt.figure()
from agent.grid import Grid, Rectangle

def CollisionDetection(agent_specs, lattice_size, dx):
	grid = Grid([lattice_size, lattice_size], dx)

	for agent, specs in agent_specs.iteritems():
		box = specs['indices']
		grid.impress(box)

	total_overlap = grid.overlap()

	return grid, total_overlap

def accept(delta):
	prob_accept = np.exp(-delta/10)

	return prob_accept


def volume_exclusion(agents):
	grid, overlap = CollisionDetection(agents, edge_length, resolution)
	while overlap > 0:
		agents_new = agents.copy()

		# update one agent at a time
		for agent_id, specs in agents_new.iteritems():
			agent = agents_new[agent_id]

			searching = True
			while searching:
				location = agent['location'] + np.random.normal(scale=np.sqrt(TRANSLATIONAL_JITTER), size=2)
				orientation = agent['orientation'] + np.random.normal(scale=ROTATIONAL_JITTER) % (2 * np.pi)

				# check if agent is in the environmental bounds
				box = Rectangle([agent['radius'], agent['length']], location, orientation)
				indices = box.render(resolution)
				searching = not grid.check_in_bounds(indices)

			agent['location'] = location
			agent['orientation'] = orientation
			agent['indices'] = indices

			grid_new, overlap_new = CollisionDetection(agents_new, edge_length, resolution)

			if overlap_new <= overlap:
				agents = agents_new
				grid = grid_new
				overlap = overlap_new

				plt.imshow(grid.grid)
				plt.pause(0.0001)

			elif np.random.rand() < accept(overlap_new - overlap):
				agents = agents_new
				grid = grid_new
				overlap = overlap_new

				plt.imshow(grid.grid)
				plt.pause(0.0001)


ROTATIONAL_JITTER = 0.1 # (radians/s)
TRANSLATIONAL_JITTER = 0.01 # (micrometers/s)

edge_length = 10
resolution = 0.1
agents = {
	'aardvark': {
		'location': (5, 5),
		'orientation': 3*np.pi/4,
		'length': 2.5,
		'radius': 0.5},
	'basilisk': {
		'location': (4, 7),
		'orientation': 0,
		'length': 4.0,
		'radius': 1.0},
	'capybara': {
		'location': (7, 6),
		'orientation': np.pi/4,
		'length': 3.0,
		'radius': 0.5},
	'dingo': {
		'location': (4, 4),
		'orientation': np.pi/5,
		'length': 4.0,
		'radius': 1.0}}

#give agents indices
for agent_id, agent in agents.iteritems():
	box = Rectangle([agent['radius'], agent['length']], agent['location'], agent['orientation'])
	indices = box.render(resolution)
	agent['indices'] = indices

volume_exclusion(agents)


import ipdb; ipdb.set_trace()