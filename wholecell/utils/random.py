#!/usr/bin/env python

"""
random.py

Special random number generators.  Most are holdovers from the original port.

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/2/2014
"""

from __future__ import division

import numpy as np

def randCounts(randomState, counts, N):
	counts = np.array(counts)
	if counts.shape == ():
		counts = counts.reshape(1)
	if np.any(counts < 0) or counts.dtype != np.dtype(np.int):
		raise Exception, "counts must contain positive integers."
	if N < 0:
		raise Exception, "N must be positive."

	cumsumCounts = np.cumsum(counts)
	positiveSelect = True

	if N > cumsumCounts[-1]:
		raise Exception, "N must be at most the total available counts."

	if N == cumsumCounts[-1]:
		return counts
	elif N > cumsumCounts[-1] / 2:
		positiveSelect = False
		N = cumsumCounts[-1] - N

	selectedCounts = np.zeros(np.shape(counts))

	for i in xrange(N):
		idx = np.ravel(np.where(randomState.randi(cumsumCounts[-1]) + 1 <= cumsumCounts))[0]
		selectedCounts[idx] += 1
		cumsumCounts[idx:] -= 1

	if not positiveSelect:
		selectedCounts = counts - selectedCounts

	return selectedCounts

def stochasticRound(randomState, value):
	# value = np.copy(valueToRound)
	value = np.array(value)
	valueShape = value.shape
	valueRavel = np.ravel(value)
	roundUp = randomState.rand(valueRavel.size) < (valueRavel % 1)
	valueRavel[roundUp] = np.ceil(valueRavel[roundUp])
	valueRavel[~roundUp] = np.floor(valueRavel[~roundUp])
	if valueShape != () and len(valueShape) > 1:
		return np.unravel_index(valueRavel, valueShape)
	else:
		return valueRavel

def make_elongation_rates_flat(
		size,
		base,
		amplified,
		ceiling,
		variable_elongation=False):
	'''
	Create an array of rates where all values are at a base rate except for a set which
	is at another rate.

	Args:
	    size (int): size of new array of rates.
	    base (int): unadjusted value for all rates.
        amplified (array[int]): indexes of each rate to adjust.
        ceiling (int): adjusted rate for amplified indexes.

	Returns:
	    rates (array[int]): new array with base and adjusted rates.
	'''

	rates = np.full(
		size,
		base,
		dtype=np.int64)

	if variable_elongation:
		rates[amplified] = ceiling

	return rates

def make_elongation_rates(
		random,
		size,
		base,
		amplified,
		ceiling,
		time_step,
		variable_elongation=False):
	'''
	Create an array of rates where all values are at a base rate except for a set which
	is at another rate. Also performs a stochastic rounding of values after applying the
	provided time step. 

	Args:
	    random (RandomState): for generating random numbers.
	    size (int): size of new array of rates.
	    base (int): unadjusted value for all rates.
        amplified (array[int]): indexes of each rate to adjust.
        ceiling (int): adjusted rate for amplified indexes.
	    time_step (float): the current time step. 
	    variable_elongation (bool): whether to add amplified values to the array.

	Returns:
	    rates (array[int]): new array with base and adjusted rates.
	'''


	rates = make_elongation_rates_flat(size, base, amplified, ceiling, variable_elongation)

	if random:
		rates = np.min([rates, stochasticRound(random, rates * time_step)], axis=0)
	else:
		rates = np.rint(rates)

	return np.array(rates, dtype=np.int64)

def randomlySelectRows(randomState, mat, prob):
	nRndRows = randomState.stochasticRound(prob * np.shape(mat)[0])
	return randomState.randomlySelectNRows(mat, nRndRows)

def randomlySelectNRows(randomState, mat, nRndRows = np.Inf):
	# mat = np.copy(matToChooseFrom)
	rndIdxs = np.sort(randomState.randsample(np.shape(mat)[0], np.min([np.shape(mat)[0], nRndRows]), False))
	mat = mat[rndIdxs, :]
	return mat, rndIdxs
