from __future__ import absolute_import, division, print_function

import numpy as np

import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
plt.ion()
fig = plt.figure()
from agent.grid import Grid, Rectangle

def CollisionDetection(agent_specs, lattice_size, dx):

	n_sites = int(lattice_size / dx)
	grid = Grid([lattice_size, lattice_size], dx)
	# collision_grid = np.zeros((n_sites, n_sites))

	for agent, specs in agent_specs.iteritems():
		location = specs['location']
		orientation = specs['orientation']
		length = specs['length']
		radius = specs['radius']

		box = Rectangle([radius, length], location, orientation)
		grid.impress(box)

		# half_length = length/2
		# cos = np.cos(orientation)
		# sin = np.sin(orientation)

		# NW = [
		# 	location[0] - half_length * cos - radius * sin,
		# 	location[1] - half_length * sin + radius * cos]

		# NE = [
		# 	location[0] + half_length * cos - radius * sin,
		# 	location[1] + half_length * sin + radius * cos]

		# SW = [
		# 	location[0] - half_length * cos + radius * sin,
		# 	location[1] - half_length * sin - radius * cos]

		# SE = [
		# 	location[0] + half_length * cos + radius * sin,
		# 	location[1] + half_length * sin - radius * cos]

		# NW_d = [int(x / dx) for x in NW]
		# NE_d = [int(x / dx) for x in NE]
		# SW_d = [int(x / dx) for x in SW]
		# SE_d = [int(x / dx) for x in SE]

		# collision_grid[NW_d[0], NW_d[1]] += 1
		# collision_grid[NE_d[0], NE_d[1]] += 1
		# collision_grid[SW_d[0], SW_d[1]] += 1
		# collision_grid[SE_d[0], SE_d[1]] += 1

	return grid



edge_length = 10
resolution = 0.1
agents = {
	'aardvark': {
		'location': (5, 5),
		'orientation': np.pi/4,
		'length': 4.0,
		'radius': 0.5},
	'basilisk': {
		'location': (5, 5),
		'orientation': 0,
		'length': 5.0,
		'radius': 1.0},
	'capybara': {
		'location': (7, 7),
		'orientation': np.pi*2/5,
		'length': 3.0,
		'radius': 2.0},
	'dingo': {
		'location': (4, 4),
		'orientation': np.pi*8/5,
		'length': 7.0,
		'radius': 1.0}}

grid = CollisionDetection(agents, edge_length, resolution)

plt.imshow(grid.grid)
plt.savefig('grid.png')
plt.pause(0.0001)

import ipdb; ipdb.set_trace()
