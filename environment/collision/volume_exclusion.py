from __future__ import absolute_import, division, print_function

import copy
import numpy as np
import argparse

# matplotlib stuff
import matplotlib
matplotlib.use('TKAgg', warn=False)
import matplotlib.pyplot as plt

from agent.grid import Grid, Rectangle

TEMPERATURE = 20 # for acceptance function


def accept(delta, temperature):
	'''
	Accept the delta at a probability dependent on the temperature.

	Args:
	    delta (float): A positive delta for gradient descent that is acceptable based on the
	        given temperature.
	    temperature (float): A higher temperature is more likely to accept the delta. 
	'''

	probability_threshold = np.exp(-delta / temperature)
	return np.random.rand() < probability_threshold


def make_shapes(agents):
	'''
	Make a dictionary out of each agent dictionary based on the function provided under its
	`render` key.

	Args:
	    agents (dict(str, dict(str, *))): Each agent dictionary requires at least a `render` key
	        with a function that accepts the agent dictionary and returns a `Shape` object.
	'''

	return {
		agent_id: agent['render'](agent)
		for agent_id, agent in agents.iteritems()}


def volume_exclusion(grid, agents, scale=1., max_cycles=100, callback=None, jitter=False):
	'''
	Perform volume exclusion on the agents in this grid.

	This function takes a grid and a dictionary of agents with locations and orientations and
	returns those same agents with their locations and orientations updated to avoid any collisions
	between them.

	Args:
	    grid (Grid): The grid these agents will be placed on.
	    agents (dict(str, dict(str, *))): A dictionary of agent dictionaries. These agents have
	        three required keys:

	        * location: the location of this agent's midpoint.
	        * orientation: the direction this agent is facing.
	        * render: a function that accepts the agent dictionary and a dx and returns a
	            `Shape` object representing the indexes into a grid of the given dx.
	'''

	shapes = make_shapes(agents)
	overlap, forces = grid.collisions(shapes)

	cycles = 0
	while overlap > 0 and cycles < max_cycles:
		potential_agents = copy.deepcopy(agents)

		for agent_id, agent in potential_agents.iteritems():
			force = forces[agent_id]

			location_jitter = 0
			orientation_jitter = 0
			if jitter:
				location_jitter = np.random.normal(scale=np.sqrt(TRANSLATIONAL_JITTER), size=2)
				orientation_jitter = np.random.normal(scale=ROTATIONAL_JITTER) % (2 * np.pi)

			agent['location'] += force * scale + location_jitter
			agent['orientation'] += orientation_jitter

		shapes_new = make_shapes(potential_agents)
		overlap_new, forces_new = grid.collisions(shapes_new)

		delta_overlap = overlap_new - overlap
		if delta_overlap <= 0 or accept(delta_overlap, TEMPERATURE):
			agents = potential_agents
			overlap = overlap_new
			shapes = shapes_new
			forces = forces_new

			if callback:
				callback(agents, overlap, forces, grid)

		cycles += 1

	return agents


if __name__ == '__main__':
	''' Perform volume exclusion on a test example and animate in matplotlib '''

	parser = argparse.ArgumentParser(description='volume exclusion')
	parser.add_argument('--animating', default=False, action='store_true')
	args = parser.parse_args()

	animating = args.animating
	if animating:
		plt.ion()
		fig = plt.figure()

	ROTATIONAL_JITTER = 0.1 # (radians/s)
	TRANSLATIONAL_JITTER = 0.0001 # (micrometers/s)

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
		'ocelot',
		'panda',
		'quail',
		'rhinoceros',
		'shark',
		'tapir',
		'urchin',
		'vole',
		'whale',
		'xenons',
		'yak',
		'zebrafish']

	def make_shape(agent):
		return Rectangle(
			[agent['radius'],
			 agent['length']],
			agent['location'],
			agent['orientation'])

	for animal in animals:
		agents[animal] = {
			'location': np.random.random(2) * 7 + 1.5,
			'orientation': np.random.random(1)[0] * np.pi * 2,
			'length': 2,
			'radius': 0.5,
			'render': make_shape}

	grid = Grid([edge_length, edge_length], resolution)

	def null_callback(agents, overlap, forces, grid):
		pass

	def animation_callback(agents, overlap, forces, grid):
		if animating:
			plt.imshow(grid.grid)
			plt.pause(0.0001)

	callback = animation_callback if animating else null_callback
	volume_exclusion(grid, agents, callback=callback)

	import ipdb; ipdb.set_trace()
