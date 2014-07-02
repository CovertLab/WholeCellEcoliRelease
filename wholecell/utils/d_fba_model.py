#!/usr/bin/env python

"""
dFbaModel.py

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 11/8/2013
"""

import numpy as np
import copy
import functools
import wholecell.utils.linear_programming

class dFbaModel(object):

	def __init__(self, metIds = None, rxns = None, mediaEx = None, biomass = None):

		self._metIds = self._createMetIds(metIds, biomass)
		self._rxnIds = self._createRxnIds(rxns, mediaEx, biomass)

		self._S = np.zeros([len(self._metIds), len(self._rxnIds)])
		self._b = np.zeros(len(self._metIds))
		self._v_upper = np.inf * np.ones(len(self._rxnIds))
		self._v_lower = -np.inf * np.ones(len(self._rxnIds))
		self._v = np.zeros(len(self._rxnIds))
		self._c = np.zeros(len(self._rxnIds))

		self._metConstTypes = ["S"] * len(self._metIds)
		self._rxnVarTypes = ["C"] * len(self._rxnIds)

		
		self._metIdToIdx = dict((x[1], x[0]) for x in enumerate(self._metIds))
		self._metIdxToId = dict((x[0], x[1]) for x in enumerate(self._metIds))

		self._rxnIdToIdx = dict((x[1], x[0]) for x in enumerate(self._rxnIds))
		self._rxnIdxToId = dict((x[0], x[1]) for x in enumerate(self._rxnIds))

		self._metGroups = {}
		self._rxnGroups = {}

		self._scaleFactor = np.NaN

		self._recalculateSolution = False

		self._createGroups(mediaEx, biomass)
		self._createScaleFactor(biomass)
		self._populateS(rxns, biomass)
		self._populateC()
		self._populateBounds()

	def _createMetIds(self, metIds, biomass):

		return \
			copy.copy(metIds)

	def _createRxnIds(self, rxns, mediaEx, biomass):

		return \
			["rxn_" + r["id"] for r in rxns] + \
			["mediaEx_" + m["rxnId"] for m in mediaEx] + \
			["g_bio"]

	def _createGroups(self, mediaEx, biomass):
		""" Create metabolite and reaction groups that we'll need (e.g., to set up the S matrix) """

		self.metGroupNew("real", [x for x in self._metIds if x[:len("pseudo")] != "pseudo"])
		self.metGroupNew("mediaEx", [x["met"] for x in mediaEx])
		self.metGroupNew("biomass", [x["id"] for x in biomass])

		self.rxnGroupNew("real", [x for x in self._rxnIds if x[:len("rnx_")] == "rxn_"])
		self.rxnGroupNew("mediaEx", [x for x in self._rxnIds if x[:len("mediaEx_")] == "mediaEx_"])
		self.rxnGroupNew("g_bio", ["g_bio"])

		self.rxnGroupNew("lowerMutable",
						self.rxnGroup("real").ids() + 
						self.rxnGroup("mediaEx").ids())
		self.rxnGroupNew("lowerImmutable",
						self.rxnGroup("g_bio").ids())
		self.rxnGroupNew("upperMutable",
						self.rxnGroup("real").ids() + 
						self.rxnGroup("mediaEx").ids() +
						self.rxnGroup("g_bio").ids())
		# self.rxnGroupNew("upperImmutable",
		# 				self.rxnGroup("x").ids() +
		# 				self.rxnGroup("(g_bio-f_i)").ids())

	def _createScaleFactor(self, biomass):
		# We want our biomass coefficients to be nicely scaled so that they are centered (in a logarithmic sense) around 1
		# That is, if our biomass coefficients were 1000, 100, and 10, we would want them to be scaled to 10, 1, and 0.1, respectively
		# This helps the numerics of the linear solver
		# _scaleFactor is the factor by which to divide the given biomass coefficients
		# For the above example, it would be 10.
		# Because we change the biomass coefficients, we also have to change the flux bounds by _scaleFactor.
		# This occurs in solution().
		# We do not alert the user to any of this. Instead, we try to handle all of this nicely behind the scenes.
		# That is, if the user uses the interface properly, they will have no idea that scaling occurred.
		# self._scaleFactor = np.exp(np.mean(np.log(np.abs(np.array([x["coeff"] for x in biomass])))))
		self._scaleFactor = np.exp(np.mean(np.log(np.abs(np.array([x["coeff"] for x in biomass if np.abs(x["coeff"]) > 1e-9])))))
		self._scaleFactor = 1.

	def _populateS(self, rxns, biomass):
		
		# Populate real reactions
		for rxn in rxns:
			rxnIdx = self.rxnIdxs(["rxn_" + rxn["id"]])[0]

			for stoich in rxn["stoichiometry"]:
				# metIdx = self.metIdxs(["%s:%s[%s]" % (stoich["molecule"], stoich["form"], stoich["location"])])[0]
				metIdx = self.metIdxs(["%s[%s]" % (stoich["molecule"], stoich["location"])])[0]

				self._S[metIdx, rxnIdx] = stoich["coeff"]

		# Populate media exchange reactions
		self._S[self.metGroup("mediaEx").idxs(), self.rxnGroup("mediaEx").idxs()] = -1.

		# Don't populate g_bio reaction (it should be all zeros)
		self._S[self.metGroup("biomass").idxs(), self.rxnGroup("g_bio").idxs()] = np.array([b["coeff"] for b in biomass])

		# S matrix has changed, so any LP solution needs to be (re)calculated
		self._recalculateSolution = True

	def _populateC(self):
		self._c[self.rxnGroup("g_bio").idxs()] = 1

		self._recalculateSolution = True

	def _populateBounds(self):
		self.v_lowerIs(self.rxnGroup("g_bio").idxs(), 0, True)

		self._recalculateSolution = True

	def metIdxs(self, metIds):
		return np.array([self._metIdToIdx[x] for x in metIds])

	def rxnIdxs(self, rxnIds):
		return np.array([self._rxnIdToIdx[x] for x in rxnIds])

	def metIds(self, metIdxs = None):
		if metIdxs == None:
			return self._metIds[:]

		return [self._metIdxToId[x] for x in metIdxs]

	def rxnIds(self, rxnIdxs = None):
		if rxnIdxs == None:
			return self._rxnIds[:]

		return [self._rxnIdxToId[x] for x in rxnIdxs]

	def metGroupNew(self, name, metIds):
		if name in self._metGroups:
			raise Exception, "Metabolite group with name '%s' already exists." % (name)

		self._metGroups[name] = metGroup(self, name, metIds)
		return self._metGroups[name]

	def metGroupDel(self, name):
		if name not in self._metGroups:
			raise Exception, "Metabolite group with name '%s' does not exist." % name

		del self._metGroups[name]

	def metGroup(self, name):
		if name not in self._metGroups:
			raise Exception, "Metabolite group with name '%s' does not exist." % name

		return self._metGroups[name]

	def metGroupNames(self):
		return [x for x in self._metGroups]

	def rxnGroupNew(self, name, rxnIds):
		if name in self._rxnGroups:
			raise Exception, "Reaction group with name '%s' already exists." % (name)

		self._rxnGroups[name] = rxnGroup(self, name, rxnIds)
		return self._rxnGroups[name]

	def rxnGroupDel(self, name):
		if name not in self._rxnGroups:
			raise Exception, "Reaction group with name '%s' does not exist." % name

		del self._rxnGroups[name]

	def rxnGroup(self, name):
		if name not in self._rxnGroups:
			raise Exception, "Reaction group with name '%s' does not exist." % name

		return self._rxnGroups[name]

	def rxnGroupNames(self):
		return [x for x in self._rxnGroups]

	def v_upper(self, idxs = None):
		if idxs == None:
			return self._v_upper

		return self._v_upper[idxs]

	def v_upperIs(self, idxs = None, values = None, overrideImmutable = False):
		if values == None:
			return

		if idxs == None:
			self._v_upper = values
		else:
			self._v_upper[idxs] = values

		# if not overrideImmutable and not np.all(self._v_upper[self.rxnGroup("upperImmutable").idxs()] == 0):
		# 	import ipdb; ipdb.set_trace()
		# 	raise Exception, "Attempting to change immutable bound."

		self._recalculateSolution = True

	def v_lower(self, idxs = None):
		if idxs == None:
			return self._v_lower

		return self._v_lower[idxs]

	def v_lowerIs(self, idxs = None, values = None, overrideImmutable = False):
		if values == None:
			return

		if idxs == None:
			self._v_lower = values
		else:
			self._v_lower[idxs] = values

		if not overrideImmutable and not np.all(self._v_lower[self.rxnGroup("lowerImmutable").idxs()] == 0):
			raise Exception, "Attempting to change immutable bound."

		self._recalculateSolution = True

	def solution(self):
		if self._recalculateSolution:

			self._scaleFactor = 1.

			v_lower = np.maximum(self._v_lower / self._scaleFactor, -1000)
			v_upper = np.minimum(self._v_upper / self._scaleFactor,  1000)

			# self._S[self.metGroup("biomass").idxs(), self.rxnGroup("g_bio").idxs()] /= self._scaleFactor

			self._v, tmp = wholecell.utils.linear_programming.linearProgramming(
			"maximize", self._c, self._S, self._b,
			v_lower, v_upper, self._metConstTypes, self._rxnVarTypes,
			None	# Will use glpk
			)

			# self._S[self.metGroup("biomass").idxs(), self.rxnGroup("g_bio").idxs()] *= self._scaleFactor

			self._recalculateSolution = False

		return self._v

	def metaboliteProduction(self, metIdxs, solution):
		pass


class metGroup(object):

	def __init__(self, model, name, metIds):
		self._model = model
		self._ids = metIds
		self._name = name

		self._idxs = model.metIdxs(metIds)

	def model(self):
		return self._model

	def ids(self):
		return self._ids

	def name(self):
		return self._name

	def idxs(self):
		return self._idxs


class rxnGroup(object):

	def __init__(self, model, name, rxnIds):
		self._model = model
		self._ids = rxnIds
		self._name = name

		self._idxs = model.rxnIdxs(rxnIds)

	def model(self):
		return self._model

	def ids(self):
		return self._ids

	def name(self):
		return self._name

	def idxs(self):
		return self._idxs



class bounds(object):

	def __init__(self, types, rxnIds, upperBound = True):
		self._upperBound = upperBound
		self._defaultValue = np.Inf
		self._mergeFunction = functools.partial(np.nanmin, axis = 0)
		if not self._upperBound:
			self._defaultValue *= -1
			self._mergeFunction = functools.partial(np.nanmax, axis = 0)

		self._bounds = self._defaultValue * np.ones((len(types), len(rxnIds)))
		self._merged = self._defaultValue * np.ones(len(rxnIds))

		self._typeToIdx = dict((x[1], x[0]) for x in enumerate(types))
		self._idxToType = dict((x[0], x[1]) for x in enumerate(types))
		
		self._rxnToIdx = dict((x[1], x[0]) for x in enumerate(rxnIds))
		self._idxToRxn = dict((x[0], x[1]) for x in enumerate(rxnIds))

		self._recalculateSolution = True

	def mergeFunction(self):
		return self._mergeFunction

	def mergeFunctionIs(self, f):
		self._mergeFunction = f

	def mergedValues(self, idxs = None):
		if self._recalculateSolution:
			self._merged = self._mergeFunction(self._bounds)
			self._recalculateSolution = False

		if idxs == None:
			return self._merged

		return self._merged[idxs]

	def values(self, idxs = None, type_ = None):
		if type_ == None:
			return

		if idxs == None:
			return self._bounds[self._typeToIdx[type_], :]

		return self._bounds[self._typeToIdx[type_], idxs]

	def valuesIs(self, idxs = None, type_ = None, values = None):
		if type_ == None or values == None:
				return

		if idxs == None:
			self._bounds[self._typeToIdx[type_], :] = values
		else:
			self._bounds[self._typeToIdx[type_], idxs] = values

		self._recalculateSolution = True