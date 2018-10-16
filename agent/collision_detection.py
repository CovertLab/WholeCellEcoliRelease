from __future__ import absolute_import, division, print_function

import copy
import numpy as np

# matplotlib stuff
import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
plt.ion()
fig = plt.figure()
from agent.grid import Grid, Rectangle

def collision_detection(grid, agents):
	grid.reset()

	shapes = {}
	forces = {}

	for agent_id, agent in agents.iteritems():
		shapes[agent_id] = agent['render'](agent)
		grid.impress(shapes[agent_id])

	overlap = grid.overlap()

	# get forces
	for agent_id, agent in agents.iteritems():
		location = agent['location']
		shape = shapes[agent_id]
		forces[agent_id] = grid.forces(location, shape)

	return overlap, forces


def accept(delta, temp):
	probability_threshold = np.exp(-delta/temp)
	return np.random.rand() < probability_threshold


def volume_exclusion(grid, agents):
	overlap, forces = collision_detection(grid, agents)
	while overlap > 0:
		potential_agents = copy.deepcopy(agents)

		# update one agent at a time
		for agent_id, agent in potential_agents.iteritems():
			force = forces[agent_id]
			location_jitter = 0 # np.random.normal(scale=np.sqrt(TRANSLATIONAL_JITTER), size=2)
			orientation_jitter = 0 # np.random.normal(scale=ROTATIONAL_JITTER) % (2 * np.pi)

			agent['location'] += force + location_jitter
			agent['orientation'] += orientation_jitter

		overlap_new, forces_new = collision_detection(grid, potential_agents)

		delta_overlap = overlap_new - overlap
		if delta_overlap <= 0 or accept(delta_overlap, TEMPERATURE):
			agents = potential_agents
			overlap = overlap_new
			forces = forces_new

			plt.imshow(grid.grid)
			plt.pause(0.0001)


ROTATIONAL_JITTER = 0.1 # (radians/s)
TRANSLATIONAL_JITTER = 0.0001 # (micrometers/s)
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

volume_exclusion(grid, agents)

import ipdb; ipdb.set_trace()
