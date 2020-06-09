#!/usr/bin/python

'''
See README.md or description.pdf.
'''

from __future__ import absolute_import, division, print_function

import os

from argparse import ArgumentParser

import numpy as np

import matplotlib.pyplot as plt

from . import parsimony
from six.moves import range

# TODO: replace some of the uglier structures with collections.namedtuple classes

# Number of model parameters - currently fixed at 7 by the model structure
# TODO: dynamically resizable problem
N_PARAMETERS = 7

# Optimization objective constants
TARGET_OUTPUT_FLUX = 1.0 # Target overall flux through pathway
TARGET_C1 = 1.0 # Target concentration of chemical species #1 (first in pathway)
MISFIT_WEIGHT = 1e-8 # Weight on target concentration misfit penalty, works best if small

# Matrix that transforms from basic (log-scaled) parameters to functional
# parameter space, in terms of (log) forward and reverse rates of reaction
TRANSFORMATION_MATRIX = np.array([
	[+1,  0,  0,  0, +1,  0,  0],
	[ 0, +1,  0,  0, +1,  0,  0],
	[ 0, +1,  0,  0,  0, +1,  0],
	[ 0,  0, +1,  0,  0, +1,  0],
	[ 0,  0, +1,  0,  0,  0, +1],
	[ 0,  0,  0, +1,  0,  0, +1],
	], np.float64)

# Search approaches
NAIVE_SEARCH_DIRECTIONS = np.identity(N_PARAMETERS)
RANDOM_SEARCH_DIRECTION_GENERATOR = lambda: random_direction_vectors(
	N_PARAMETERS, N_PARAMETERS
	)
PROTOPARSIMONIOUS_SEARCH_DIRECTIONS = np.linalg.pinv(
	TRANSFORMATION_MATRIX
	).transpose()
PARSIMONIOUS_SEARCH_DIRECTIONS = parsimony.construct_parsimony_objects(
	np.vstack([TRANSFORMATION_MATRIX, np.identity(N_PARAMETERS)]),
	TRANSFORMATION_MATRIX
	)[0]

# Iterable over the approaches - comment out any line you're not interested in
SEARCH_APPROACHES = [
	('Naive', NAIVE_SEARCH_DIRECTIONS),
	('Random', RANDOM_SEARCH_DIRECTION_GENERATOR),
	('Proto-parsimonious', PROTOPARSIMONIOUS_SEARCH_DIRECTIONS),
	('Parsimonious', PARSIMONIOUS_SEARCH_DIRECTIONS),
	]

N_APPROACHES = len(SEARCH_APPROACHES)

# Initial parameter value scale as well as optimization sweep scale, in base-10
SCALE = 2


# Optimization algorithm constants

MAX_ITERATIONS = 2000
SAMPLES = 31 # Number of samples to take in each direction - should be odd

# Scale adjustment parameters

# Weight on new scale vs. old scale
# Larger = more aggressive
WEIGHT_NEW_SCALE = 1.0

# Bias towards larger scales - important if current search scale is too small
# Scale will stall out or explode if too large
UPWARDS_BIAS = 0.1

# Below these values, the algorithm will stop (to save iterations)
MIN_OBJECTIVE_VALUE = 1e-15
# Possible improvement: quit if objective (fold) change is too small for too long
MIN_SCALE = 1e-20 # usually not triggered


# Plotting

# Threshold for declaring a successful optimization
FIGSIZE = (5, 5)
TARGET_OBJECTIVE_VALUE = 1e-2 * MISFIT_WEIGHT

OUTDIR = os.path.join(
	os.path.split(__file__)[0],
	'out'
	)

# Functions

# TODO: try to reduce the excessive transposition logic

def model(parameter_values):
	'''
	This is very explicit to give a clear idea of what this model is.
	'''

	c1, c2, c3, c4, kA, kB, kC = 10**(parameter_values)

	vA = kA * (c1 - c2)
	vB = kB * (c2 - c3)
	vC = kC * (c3 - c4)

	v = np.array([vA, vB, vC]).transpose()

	dc2 = vB - vA
	dc3 = vC - vB

	dc_dt = np.array([dc2, dc3]).transpose()

	dln_c2 = dc2 / c2
	dln_c3 = dc3 / c3

	dln_c_dt = np.array([dln_c2, dln_c3]).transpose()

	return (v, dc_dt, dln_c_dt)

def objective_function(parameter_values):
	(v, dc_dt, dln_c_dt) = model(parameter_values)

	absolute_deviation_from_ss = np.sum(np.square(dc_dt.transpose()), 0)
	relative_deviation_from_ss = np.sum(np.square(dln_c_dt.transpose()), 0)

	misfit = np.abs(parameter_values[0] - np.log10(TARGET_C1))

	target_output_flux_deviation = np.square(
		v.transpose()[-1] - TARGET_OUTPUT_FLUX
		)

	return (
		absolute_deviation_from_ss
		+ relative_deviation_from_ss
		+ target_output_flux_deviation
		+ MISFIT_WEIGHT*misfit
		)

def optimization_step(
		initial_parameter_values,
		objective_function,
		scale,
		number_of_samples,
		search_directions
		):

	reference_parameter_stack = np.tile(
		initial_parameter_values, (number_of_samples, 1)
		).transpose()

	shifts = np.linspace(-scale, +scale, number_of_samples)

	best_direction = None
	best_shift = None
	best_objective_value = objective_function(initial_parameter_values)

	for (search_direction_index, search_direction) in enumerate(search_directions):
		parameter_stack = reference_parameter_stack.copy()

		parameter_stack += search_direction[:, None] * shifts[None, :]

		perturbed_objective_values = objective_function(parameter_stack)

		# Possible improvement: bounding
		# This naive bounding approach tends to get 'stuck' pretty badly
		# perturbed_objective_values[np.any(parameter_stack < -SCALE, 0)] = np.nan
		# perturbed_objective_values[np.any(parameter_stack > +SCALE, 0)] = np.nan

		smallest_where = np.nanargmin(perturbed_objective_values)
		smallest_objective_value = perturbed_objective_values[smallest_where]

		if smallest_objective_value < best_objective_value:
			best_direction = search_direction_index
			best_shift = shifts[smallest_where]
			best_objective_value = smallest_objective_value

	best_parameters = initial_parameter_values.copy()

	if best_shift is None:
		best_shift = 0.0

	else:
		best_parameters += best_shift * search_directions[best_direction]

	return (best_parameters, best_objective_value, best_shift, best_direction)

def random_direction_vectors(number, size):
	vectors = np.random.normal(size = (number, size))
	vectors /= np.sqrt(np.sum(np.square(vectors), 1))[:, None]

	return vectors

def call_if_callable(obj):
	if callable(obj):
		out = obj()

	else:
		out = obj

	return out

def main():
	# Initialize
	argparser = ArgumentParser()
	argparser.add_argument(
		'-s', '--seed',
		type = int, default = None,
		help = (
			'The seed used for problem initialization (as well as for '
			+ 'generating the random search vectors).  If not provided, a '
			+ 'seed will be selected randomly.')
		)
	argparser.add_argument(
		'-i', '--interactive',
		action = 'store_true',
		help = (
			'Enables interactive mode.  Otherwise the plots are saved to the '
			+'out/ directory.'
			)
		)

	args = argparser.parse_args()

	if args.seed is None:
		seed = np.random.randint(2**32-1)

	else:
		seed = args.seed

	print('Initializing with random seed {}'.format(seed))
	# TODO: pass around a proper numpy.random.RandomStream object
	np.random.seed(seed)

	# TODO: uniform initialization option
	initial_parameter_values = np.random.uniform(
		-SCALE, +SCALE,
		size = N_PARAMETERS
		)
	initial_objective_value = objective_function(initial_parameter_values)

	(fig_objective, axes_objective) = plt.subplots(figsize = FIGSIZE)
	(fig_scale, axes_scale) = plt.subplots(figsize = FIGSIZE)
	(fig_directions, axes_directions) = plt.subplots(figsize = FIGSIZE, nrows = N_APPROACHES)
	(fig_direction_vectors, axes_direction_vectors) = plt.subplots(figsize = FIGSIZE, nrows = N_APPROACHES)

	if N_APPROACHES == 1:
		axes_directions = (axes_directions,)

	# Optimize using each approach
	for (approach_index, (name, search_directions)) in enumerate(SEARCH_APPROACHES):

		# TODO: functionalize inner loop
		parameter_values = initial_parameter_values.copy()
		objective_values = [initial_objective_value]
		scale = SCALE
		scales = [scale]
		directions = np.zeros((MAX_ITERATIONS, N_PARAMETERS), np.float64)
		direction_vectors = np.zeros(
			(MAX_ITERATIONS, call_if_callable(search_directions).shape[0]),
			np.bool
			)

		for iteration in range(1, MAX_ITERATIONS+1):
			sd = call_if_callable(search_directions)

			(
				new_parameter_values,
				new_objective_value,
				best_shift,
				best_direction
				) = optimization_step(
					parameter_values,
					objective_function,
					scale,
					SAMPLES,
					sd
					)

			objective_values.append(new_objective_value)
			parameter_values = new_parameter_values

			scale = (
				(scale + WEIGHT_NEW_SCALE * np.abs(best_shift))
				/ (1 + WEIGHT_NEW_SCALE)
				)
			scale *= 1 + UPWARDS_BIAS

			# Possible improvement: reset scale if it becomes too small

			scales.append(scale)

			if best_direction is not None:

				directions[iteration-1] = sd[best_direction] * np.sign(best_shift)
				direction_vectors[iteration-1, best_direction] = True

			if (
					(objective_values[-1] < MIN_OBJECTIVE_VALUE)
					or (scale < MIN_SCALE)
					):
				break

		axes_objective.plot(objective_values, label = name)
		axes_scale.plot(scales, label = name)

		directions /= np.nanmax(np.abs(directions))

		ad = axes_directions[approach_index]

		ad.imshow(
			directions.transpose(),
			aspect = 'auto',
			interpolation = 'nearest',
			# These limits keep the extremes from being too dark (and thus,
			# indistinguishable)
			vmin = -1.4, vmax = +1.4,
			cmap = 'RdBu',
			)
		ad.set_ylabel(name, fontsize = 8)
		ad.set_xticks([])
		ad.set_yticks([])

		adv = axes_direction_vectors[approach_index]

		adv.imshow(
			direction_vectors.transpose(),
			aspect = 'auto',
			interpolation = 'nearest',
			cmap = 'bone_r',
			)
		adv.set_ylabel(name, fontsize = 8)
		adv.set_xticks([])
		adv.set_yticks([])

	# Finalize figures

	axes_objective.axhline(
		TARGET_OBJECTIVE_VALUE,
		c = 'k', lw = 0.5,
		label = 'Acceptable solution threshold'
		)
	axes_objective.axhline(
		MIN_OBJECTIVE_VALUE,
		c = 'k', lw = 1.0, ls = ':',
		label = 'Early termination threshold'
		)
	axes_objective.set_yscale('log')

	# axes_objective.set_title('Objective values')

	fig_objective.tight_layout()

	axes_objective.legend(loc = 'best')

	# Don't plot, since it usually doesn't come into play
	# axes_scale.axhline(
	# 	MIN_SCALE,
	# 	c = 'k', lw = 1.0, ls = ':',
	# 	label = 'Early termination threshold'
	# 	)

	axes_scale.set_yscale('log')

	# axes_scale.set_title('Search scale')

	fig_scale.tight_layout()

	axes_scale.legend(loc = 'best')

	# fig_directions.suptitle('Search directions')
	fig_directions.tight_layout()

	# fig_direction_vectors.suptitle('Search direction vectors')
	fig_direction_vectors.tight_layout()

	if args.interactive:
		plt.show()

	else:
		outdir = os.path.join(OUTDIR, 'seed{}'.format(seed))

		if not os.path.exists(outdir):
			os.makedirs(outdir)

		for (i, (figure, name)) in enumerate([
				(fig_objective, 'objective'),
				(fig_scale, 'scale'),
				(fig_directions, 'directions'),
				(fig_direction_vectors, 'direction_vectors'),
				]):

			filename = 'figure{}-{}.pdf'.format(i+1, name)

			figure.savefig(os.path.join(outdir, filename))

			plt.close(figure)

if __name__ == '__main__':
	main()
