from __future__ import absolute_import, division, print_function

import numpy as np

import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
plt.ion()
fig = plt.figure()


def CollisionDetection(agent_specs, lattice_size, dx):

	n_sites = int(lattice_size / dx)
	collision_grid = np.zeros((n_sites, n_sites))

	for agent, specs in agent_specs.iteritems():
		location = specs['location']
		orientation = specs['orientation']
		length = specs['length']
		radius = specs['radius']
		half_length = length/2

		NW = [location[0] - half_length * np.cos(orientation) - radius * np.sin(orientation),
				   location[1] - half_length * np.sin(orientation) + radius * np.cos(orientation)]

		NE = [location[0] + half_length * np.cos(orientation) - radius * np.sin(orientation),
				   location[1] + half_length * np.sin(orientation) + radius * np.cos(orientation)]

		SW = [location[0] - half_length * np.cos(orientation) + radius * np.sin(orientation),
				   location[1] - half_length * np.sin(orientation) - radius * np.cos(orientation)]

		SE = [location[0] + half_length * np.cos(orientation) + radius * np.sin(orientation),
				   location[1] + half_length * np.sin(orientation) - radius * np.cos(orientation)]

		NW_d = [int(x / dx) for x in NW]
		NE_d = [int(x / dx) for x in NE]
		SW_d = [int(x / dx) for x in SW]
		SE_d = [int(x / dx) for x in SE]

		collision_grid[NW_d[0], NW_d[1]] += 1
		collision_grid[NE_d[0], NE_d[1]] += 1
		collision_grid[SW_d[0], SW_d[1]] += 1
		collision_grid[SE_d[0], SE_d[1]] += 1

	return collision_grid



edge_length = 10
resolution = 0.1
agents = {1 : {'location' : (5,5), 'orientation' : np.pi/4, 'length' : 4.0, 'radius' : 0.5}}

grid = CollisionDetection(agents, edge_length, resolution)

plt.imshow(grid)
plt.pause(0.0001)

import ipdb; ipdb.set_trace()