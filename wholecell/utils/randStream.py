#!/usr/bin/env python

"""
RandStream
Wrapper over numpy.random.RandomState

Adds method:
 - stochasticRound: rounding of floats to integers weighted by decimal part

Note:
 - any size arguments from Matlab need to be carefully checked for intended behavior
   - "rand(10)" in Matlab returns a 10x10 matrix
   - "numpy.random.rand(10)" returns an array with 10 elements
 - modifications have been made for 0-indexing in Python

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/22/2013
"""

import numpy
import numpy.random

class RandStream(object):
	""" Random Stream """

	def __init__(self, seed = None):
		self.randStream = numpy.random.RandomState()
		self.seed = seed
		self.__initialState__ = self.randStream.get_state()
		self.defaultStream = None

	def reset(self, seed = None):
		self.randStream.set_state(self.__initialState__)
		if seed:
			self.seed = seed

	# def setDefault(self):
	# 	# Store the default (global) random state
	# 	self.defaultStream = numpy.random.get_state()

	# 	# Set the global random state to have this random state
	# 	numpy.random.set_state(self.randStream.get_state())

	# def resetDefault(self):
	# 	numpy.random.set_state(self.defaultStream)

	def rand(self, *args):
		return self.randStream.rand(*args)

	# TODO: Change name to "randint"?
	def randi(self, imax, size = None):
		return self.randStream.randint(0, imax, size)
		# return self.randStream.randint(1, imax + 1, size)

	def randn(self, *args):
		return self.randStream.randn(*args)

	def randperm(self, n):
		array = numpy.arange(n)
		# array = numpy.arange(1, n + 1)
		self.randStream.shuffle(array)
		return array

	def random(self):
		raise Exception, "No numpy equivalent. Call appropriate function in numpy.random instead."

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
		return self.randStream.multivariate_normal(mean, cov, size)

	def randsample(self, n, k, replacement = False, w = None):
		k = numpy.array(k, dtype = numpy.int)
		if k.size != 1:
			raise Exception, "Expect k to be a scalar"

		if w != None:
			# The w argument can be weights, but choice() needs strict probabilities
			w_array = numpy.array(w, dtype = numpy.float)
			if numpy.any(w_array < 0):
				raise ValueError, "probabilities are not non-negative"
			if numpy.any(w_array == numpy.Inf) or numpy.any(w_array == numpy.NaN):
				raise ValueError, "expect finite weights"
			if numpy.abs(numpy.sum(w_array)) < 1e-9:
				raise ValueError, "weights sum to zero - cannot scale"
			w_array /= numpy.sum(w_array)
		else:
			w_array = None

		return self.randStream.choice(numpy.arange(n), k, replacement, w_array)
		# return self.randStream.choice(numpy.arange(1, n + 1), k, replacement, w)

	def randCounts(self, counts, N):
		counts = numpy.array(counts)
		if counts.shape == ():
			counts = counts.reshape(1)
		if numpy.any(counts < 0) or counts.dtype != numpy.dtype(numpy.int):
			raise Exception, "counts must contain positive integers."
		if N < 0:
			raise Exception, "N must be positive."

		cumsumCounts = numpy.cumsum(counts)
		positiveSelect = True

		if N > cumsumCounts[-1]:
			raise Exception, "N must be at most the total available counts."

		if N == cumsumCounts[-1]:
			return counts
		elif N > cumsumCounts[-1] / 2:
			positiveSelect = False
			N = cumsumCounts[-1] - N

		selectedCounts = numpy.zeros(numpy.shape(counts))

		for i in xrange(N):
			idx = numpy.ravel(numpy.where(self.randi(cumsumCounts[-1]) + 1 <= cumsumCounts))[0]
			selectedCounts[idx] += 1
			cumsumCounts[idx:] -= 1

		if not positiveSelect:
			selectedCounts = counts - selectedCounts

		return selectedCounts

	def stochasticRound(self, value):
		# value = numpy.copy(valueToRound)
		value = numpy.array(value)
		valueShape = value.shape
		valueRavel = numpy.ravel(value)
		roundUp = self.rand(valueRavel.size) < (valueRavel % 1)
		valueRavel[roundUp] = numpy.ceil(valueRavel[roundUp])
		valueRavel[~roundUp] = numpy.floor(valueRavel[~roundUp])
		if valueShape != () and len(valueShape) > 1:
			return numpy.unravel_index(valueRavel, valueShape)
		else:
			return valueRavel

	def randomlySelectRows(self, mat, prob):
		nRndRows = self.stochasticRound(prob * numpy.shape(mat)[0])
		return self.randomlySelectNRows(mat, nRndRows)

	def randomlySelectNRows(self, mat, nRndRows = numpy.Inf):
		# mat = numpy.copy(matToChooseFrom)
		rndIdxs = numpy.sort(self.randsample(numpy.shape(mat)[0], numpy.min([numpy.shape(mat)[0], nRndRows]), False))
		mat = mat[rndIdxs, :]
		return mat, rndIdxs

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
			numpy.array_equal(stateSelf[1], stateOther[1]) and \
			numpy.array_equal(stateSelf[2], stateOther[2]) and \
			numpy.array_equal(stateSelf[3], stateOther[3])