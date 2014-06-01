#!/usr/bin/env python

"""
RandStream
Wrapper over np.random.RandomState

Adds method:
 - stochasticRound: rounding of floats to integers weighted by decimal part

Note:
 - any size arguments from Matlab need to be carefully checked for intended behavior
   - "rand(10)" in Matlab returns a 10x10 matrix
   - "np.random.rand(10)" returns an array with 10 elements
 - modifications have been made for 0-indexing in Python

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/22/2013
"""

from __future__ import division

import numpy as np
import warnings

class RandStream(object):
	""" Random Stream """

	def __init__(self, seed = None):
		self.randStream = np.random.RandomState()
		self.seed = seed
		self.__initialState__ = self.randStream.get_state()
		self.defaultStream = None

	def reset(self, seed = None):
		self.randStream.set_state(self.__initialState__)
		if seed:
			self.seed = seed

	# def setDefault(self):
	# 	# Store the default (global) random state
	# 	self.defaultStream = np.random.get_state()

	# 	# Set the global random state to have this random state
	# 	np.random.set_state(self.randStream.get_state())

	# def resetDefault(self):
	# 	np.random.set_state(self.defaultStream)

	def rand(self, *args):
		return self.randStream.rand(*args)

	# TODO: Change name to "randint"?
	def randi(self, imax, size = None):
		return self.randStream.randint(0, imax, size)
		# return self.randStream.randint(1, imax + 1, size)

	def randn(self, *args):
		return self.randStream.randn(*args)

	def randperm(self, n):
		array = np.arange(n)
		# array = np.arange(1, n + 1)
		self.randStream.shuffle(array)
		return array

	def random(self):
		raise Exception, "No np equivalent. Call appropriate function in np.random instead."

	# TODO: Change name to "poisson"?
	def poissrnd(self, Lambda, size = None):
		return self.randStream.poisson(Lambda, size)

	# TODO: Change name to "multinomial"?
	def mnrnd(self, n, pvals, size = None):
		return self.randStream.multinomial(n, pvals, size)

	# TODO: Change name to "binomial"?
	def binornd(self, n, p, size = None):
		return self.randStream.binomial(n, p, size)

	def normal(self, loc = 0.0, scale = 1.0 , size = None):
		return self.randStream.normal(loc, scale, size)

	def multivariate_normal(self, mean, cov, size = None):
		with warnings.catch_warnings():
			warnings.filterwarnings("ignore", category=DeprecationWarning)
			value = self.randStream.multivariate_normal(mean, cov, size)

		return value

	def randsample(self, n, k, replacement = False, w = None):
		k = np.array(k, dtype = np.int)
		if k.size != 1:
			raise Exception, "Expect k to be a scalar"

		if w != None:
			# The w argument can be weights, but choice() needs strict probabilities
			w_array = np.array(w, dtype = np.float)
			if np.any(w_array < 0):
				raise ValueError, "probabilities are not non-negative"
			if np.any(w_array == np.Inf) or np.any(w_array == np.NaN):
				raise ValueError, "expect finite weights"
			if np.abs(np.sum(w_array)) < 1e-9:
				raise ValueError, "weights sum to zero - cannot scale"
			w_array /= np.sum(w_array)
		else:
			w_array = None

		return self.randStream.choice(np.arange(n), k, replacement, w_array)
		# return self.randStream.choice(np.arange(1, n + 1), k, replacement, w)

	def randCounts(self, counts, N):
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
			idx = np.ravel(np.where(self.randi(cumsumCounts[-1]) + 1 <= cumsumCounts))[0]
			selectedCounts[idx] += 1
			cumsumCounts[idx:] -= 1

		if not positiveSelect:
			selectedCounts = counts - selectedCounts

		return selectedCounts

	def stochasticRound(self, value):
		# value = np.copy(valueToRound)
		value = np.array(value)
		valueShape = value.shape
		valueRavel = np.ravel(value)
		roundUp = self.rand(valueRavel.size) < (valueRavel % 1)
		valueRavel[roundUp] = np.ceil(valueRavel[roundUp])
		valueRavel[~roundUp] = np.floor(valueRavel[~roundUp])
		if valueShape != () and len(valueShape) > 1:
			return np.unravel_index(valueRavel, valueShape)
		else:
			return valueRavel

	def randomlySelectRows(self, mat, prob):
		nRndRows = self.stochasticRound(prob * np.shape(mat)[0])
		return self.randomlySelectNRows(mat, nRndRows)

	def randomlySelectNRows(self, mat, nRndRows = np.Inf):
		# mat = np.copy(matToChooseFrom)
		rndIdxs = np.sort(self.randsample(np.shape(mat)[0], np.min([np.shape(mat)[0], nRndRows]), False))
		mat = mat[rndIdxs, :]
		return mat, rndIdxs

	def numpyShuffle(self, *args, **kwargs):
		self.randStream.shuffle(*args, **kwargs)

	def numpyChoice(self, *args, **kwargs):
		self.randStream.choice(*args, **kwargs)

	@property
	def type(self):
		return self.randStream.get_state()[0]

	@property
	def seed(self):
		pass

	@seed.setter
	def seed(self, value):
		self.randStream.seed(value)

	@property
	def state(self):
		return self.randStream.get_state()

	@state.setter
	def state(self, value):
		self.randStream.set_state(value)

	def __eq__(self, other):
		stateSelf = self.state
		stateOther = other.state

		return (stateSelf[0] == stateOther[0]) and \
			np.array_equal(stateSelf[1], stateOther[1]) and \
			np.array_equal(stateSelf[2], stateOther[2]) and \
			np.array_equal(stateSelf[3], stateOther[3])