#!/usr/bin/python

'''

A demonstration of a simple parameter estimation procedure on a standard
Michaelis-Menten rate law, with a desired rate of reaction and no prior
observed parameter values.  The system has been reduced and transformed into
two parameters for better visualization (but without loss of generality).

TODO: shrink color bars, better ticks

'''

from __future__ import absolute_import, division, print_function

import os
from typing import Any

import numpy as np

import matplotlib.pyplot as plt


# Constants (independent and derived)

# Optimization constants

TARGET_RATE_OF_REACTION = 1.0

# Parameter sweep constants

# The parameter search window (up and down).  Should be small and positive, no
# less than 1e-15 (float64 precision limit) and no greater than 1.
SPREAD = 1e-3

ALPHA_MIN = np.log(TARGET_RATE_OF_REACTION * SPREAD)
ALPHA_MAX = np.log(TARGET_RATE_OF_REACTION * SPREAD**-1)

BETA_MIN = np.log(SPREAD)
BETA_MAX = np.log(SPREAD**-1)

N_POINTS = 1000 # Needs to be large to see minima

# Plotting constants
IMSHOW_STANDARD_KWARGS = dict(
	interpolation = 'nearest', # "nearest" = no interpolation (nearest neighbor)
	origin = 'lower', # place origin in bottom-left (standard for plotting)
	)

# Color maps
CMAP_GENERAL = 'viridis' # arbitrary values, "zero" has no particular meaning
CMAP_DIVERGENT = 'RdBu' # positive/negative values, "zero" has a specific meaning (must be centered!)
CMAP_ONESIDED = 'bone' # minimum value has a specific meaning
CMAP_ROTATIONAL = 'hsv' # circularly permuted values (min/max must be chosen appropriately)

EXTENT = (
	ALPHA_MIN, ALPHA_MAX,
	BETA_MIN, BETA_MAX
	)

N_TICKS = 5 # should be a small integer, preferably odd

TICKS_FRACTIONAL_MARGIN = 0.9 # prevent ticks from being too close to the edge

FIGSIZE = (4, 4)

OUTDIR = os.path.join(
	os.path.split(__file__)[0],
	'out'
	)
EXTENSION = ( # Choose one
	# 'png'
	'pdf'
	)

# 'Raw' string needed to avoid clashes between escape characters for Python and
# for the LaTeX renderer.
ALPHA_STRING = r'$\alpha$'
BETA_STRING = r'$\beta$'

# Functions

# Model functions

def rate_of_reaction(alpha, beta):
	'''
	Compute the rate of reaction subject to two lumped parameters.

	Parameters:
		alpha: float
			Natural logarithm of the maximum velocity (AKA rate of reaction).

		beta: float
			Natural logarithm of the saturation ratio (metabolite concetration
			divided by corresponding Michaelis-Menten saturation constant).

	Returns:
		The "velocity" or rate of reaction. (float)

	Formulation:

	Starting from standard MM kinetics,

	v = k E C / (C + K)

	v: rate of reaction
	k: catalytic rate constant
	E: enzyme concentration
	C: metabolite concentration
	K: Michaelis-Menten saturation constant

	Lump the parameters:

	vmax = k E: maximum "velocity" or rate of reaction
	s = C/K: saturation ratio

	Now, rewrite the system,

	v = vmax * s / (s + 1)

	Without loss of generality, rewrite in terms of logarithmic parameters:

	alpha = ln(vmax)
	beta = ln(s)

	Now, finally,

	v = exp(alpha + beta) / (exp(beta) + 1)

	Finally, for numerical purposes, we rewrite as

	v = exp(alpha) / (1 + exp(-beta))

	Apart from their physical interpretation, alpha can be thought of as the
	height of the curve, and beta can be thought of as the horizontal shift.
	'''
	return np.exp(alpha) / (1 + np.exp(-beta))


# Optimization functions

def objective(alpha, beta):
	'''
	Returns the square of the difference between the target rate of reaction
	and the calculated rate of reaction based on the two parameters alpha and
	beta.  See the documentation for 'rate_of_reaction' for more information.
	'''
	return (rate_of_reaction(alpha, beta) - TARGET_RATE_OF_REACTION)**2

def grad(alpha, beta):
	'''
	Computes the gradient of the objective with respect to the model
	parameters.

	Note that these derivatives were computed by hand.  As such I'm not
	perfectly confident in them.
	'''
	# TODO: validate

	ror = rate_of_reaction(alpha, beta)
	dobj_dalpha = 2 * ror * (ror - TARGET_RATE_OF_REACTION)
	dobj_dbeta = dobj_dalpha / (np.exp(beta) + 1)

	return dobj_dalpha, dobj_dbeta


# Plotting functions

def imshow_percentile(array, bottom, top, **kwargs):
	'''
	Plot an array as an image, choosing upper and lower limits using
	pecentiles.
	'''
	q = [bottom, top]  # type: Any
	(vmin, vmax) = np.percentile(array, q)

	plt.imshow(array, vmin = vmin, vmax = vmax, **kwargs)

def imshow_percentile_centered(array, top, **kwargs):
	'''
	Plot an array as an image, choosing upper and lower limits using
	a percentile.  Data will be centered about zero.
	'''
	vmax = np.percentile(np.abs(array), top)
	vmin = -vmax

	plt.imshow(array, vmin = vmin, vmax = vmax, **kwargs)

def choose_ticks(min_value, max_value, n_ticks, fractional_margin):
	'''
	Select plotted ticks.
	'''
	return np.linspace(
		np.ceil(min_value * fractional_margin),
		np.floor(max_value * fractional_margin),
		n_ticks
		)

def ax(ticks, label, tick_function, label_function, ticks_only = False):
	'''
	Axis formatting function.
	'''
	if ticks_only:
		tick_function(ticks, ['']*ticks.size)

	else:
		tick_function(ticks)
		label_function(label)

def finalize(figure_number, title):
	plt.tight_layout()
	plt.savefig(os.path.join(OUTDIR,
		'figure{}-{}.{}'.format(figure_number, title, EXTENSION)
		))
	plt.clf()

# Execution

def main():
	# Sweep the parameter space, computing the objective and gradient

	alpha_range = np.linspace(ALPHA_MIN, ALPHA_MAX, N_POINTS)
	beta_range = np.linspace(BETA_MIN, BETA_MAX, N_POINTS)

	(alpha, beta) = np.meshgrid(alpha_range, beta_range)

	obj = objective(alpha, beta)

	(grad_alpha, grad_beta) = grad(alpha, beta)

	# Plotting

	if not os.path.exists(OUTDIR):
		os.mkdir(OUTDIR)

	# The derivative of a logarithm of a function is the derivative of the function
	# divided by the function itself.
	log_grad_alpha = grad_alpha / obj
	log_grad_beta = grad_beta / obj

	xticks = choose_ticks(
		ALPHA_MIN, ALPHA_MAX,
		N_TICKS, TICKS_FRACTIONAL_MARGIN
		)
	yticks = choose_ticks(
		BETA_MIN, BETA_MAX,
		N_TICKS, TICKS_FRACTIONAL_MARGIN
		)

	# Conveience functions - defined here because they refer to values that
	# strictly exist in this scope.

	def xax(ticks_only = False):
		ax(xticks, ALPHA_STRING, plt.xticks, plt.xlabel, ticks_only)

	def yax(ticks_only = False):
		ax(yticks, BETA_STRING, plt.yticks, plt.ylabel, ticks_only)


	# Figure 1: Raw objective value

	# The raw objective value is hard to see, but is provided for context.

	plt.figure(figsize = FIGSIZE)

	imshow_percentile(
		obj,
		0, 95, # Value grows very rapidly - use 95th percentile as max
		**dict(
			IMSHOW_STANDARD_KWARGS,
			cmap = CMAP_ONESIDED,
			extent = EXTENT
			)
		)

	plt.colorbar()

	xax()
	yax()

	finalize(1, 'raw_objective_values')


	# Figure 2: Log-scaled objective value

	# Log-scaling the objective value makes the features of the objective value
	# space clear.

	plt.figure(figsize = FIGSIZE)

	imshow_percentile(
		np.log(obj),
		1, 100, # objective falls off to zero at minima
		**dict(
			IMSHOW_STANDARD_KWARGS,
			cmap = CMAP_ONESIDED,
			extent = EXTENT
			)
		)

	plt.colorbar()

	xax()
	yax()

	finalize(2, 'log_objective_values')


	# Figure 3: Derivative of log objective w.r.t alpha

	plt.figure(figsize = FIGSIZE)

	imshow_percentile_centered(
		log_grad_alpha,
		95, # log gradient near minima is very steep since the objective approaches zero
		**dict(
			IMSHOW_STANDARD_KWARGS,
			cmap = CMAP_DIVERGENT,
			extent = EXTENT
			)
		)

	plt.colorbar()

	xax()
	yax()

	finalize(3, 'derivative_alpha')


	# Figure 4: Derivative of log objective w.r.t beta

	plt.figure(figsize = FIGSIZE)

	imshow_percentile_centered(
		log_grad_beta,
		95, # log gradient near minima is very steep since the objective approaches zero
		**dict(
			IMSHOW_STANDARD_KWARGS,
			cmap = CMAP_DIVERGENT,
			extent = EXTENT
			)
		)

	plt.colorbar()

	xax()
	yax()

	finalize(4, 'derivative_beta')


	# Figure 5: Optimal search direction

	plt.figure(figsize = FIGSIZE)

	# Note that the best search angle doesn't change if we are minimizing the
	# objective or the logarithmic objective, since they differ only by a scaling
	# factor.
	best_search_angle = np.arctan2(-grad_beta, -grad_alpha) # Negative because we are minimizing

	plt.imshow(
		best_search_angle,
		vmin = -np.pi, vmax = +np.pi, # Range of arctan2
		cmap = CMAP_ROTATIONAL,
		extent = EXTENT,
		**IMSHOW_STANDARD_KWARGS
		)

	plt.colorbar()

	xax()
	yax()

	finalize(5, 'best_local_search_direction')


	# Figure 6: Flow plot

	plt.figure(figsize = FIGSIZE)

	plt.streamplot(alpha, beta, -grad_alpha, -grad_beta, density = [0.5, 0.5])
	plt.axis('square')

	xax()
	yax()

	# Need to set limits explicitly, unlike with plt.imshow
	plt.xlim(ALPHA_MIN, ALPHA_MAX)
	plt.ylim(BETA_MIN, BETA_MAX)

	finalize(6, 'gradient_descent_flow')

if __name__ == '__main__':
	main()
