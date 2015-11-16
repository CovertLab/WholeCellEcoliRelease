#!/usr/bin/env python

"""
FlexTFbaModel.py

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 10/15/2013
"""

from __future__ import division

import numpy as np
import copy
import functools
import wholecell.utils.linear_programming

class FlexTFbaModel(object):

	def __init__(self, metIds = None, rxns = None, mediaEx = None, biomass = None, atpId = "ATP:mature[c]", params = None):

		self._metIds = self._createMetIds(metIds, biomass, atpId)
		self._rxnIds = self._createRxnIds(rxns, mediaEx, biomass, atpId)

		self._S = np.zeros([len(self._metIds), len(self._rxnIds)])
		self._b = np.zeros(len(self._metIds))
		self._v_upper = np.inf * np.ones(len(self._rxnIds))
		self._v_lower = -np.inf * np.ones(len(self._rxnIds))
		self._v = np.zeros(len(self._rxnIds))
		self._c = np.zeros(len(self._rxnIds))

		self._metConstTypes = ["S"] * len(self._metIds)
		self._rxnVarTypes = ["C"] * len(self._rxnIds)

		try:
			self._alpha = params["alpha"]
			self._beta = params["beta"]
			self._gamma = params["gamma"]
		except:
			self._alpha = 10.
			self._beta = 1000.
			self._gamma = -1.
		
		self._metIdToIdx = dict((x[1], x[0]) for x in enumerate(self._metIds))
		self._metIdxToId = dict((x[0], x[1]) for x in enumerate(self._metIds))

		self._rxnIdToIdx = dict((x[1], x[0]) for x in enumerate(self._rxnIds))
		self._rxnIdxToId = dict((x[0], x[1]) for x in enumerate(self._rxnIds))

		self._metGroups = {}
		self._rxnGroups = {}

		self._scaleFactor = np.NaN

		self._recalculateSolution = False

		self._createGroups(mediaEx, biomass, atpId)
		self._createScaleFactor(biomass)
		self._populateS(rxns, biomass)
		self._populateC()
		self._populateBounds()

	def _createMetIds(self, metIds, biomass, atpId):

		return \
			copy.copy(metIds) + \
			["pseudo_(f_atp-" + b["id"] + ")_" for b in biomass if b["id"] != atpId] + \
			["pseudo_(g_bio-" + b["id"] + ")_" for b in biomass]

	def _createRxnIds(self, rxns, mediaEx, biomass, atpId):

		return \
			["rxn_" + r["id"] for r in rxns] + \
			["mediaEx_" + m["rxnId"] for m in mediaEx] + \
			["f_" + b["id"] for b in biomass] + \
			["x_" + b["id"] for b in biomass] + \
			["(f_atp-f_" + b["id"] + ")" for b in biomass if b["id"] != atpId] + \
			["g_bio"] + \
			["(g_bio-f_" + b["id"] + ")" for b in biomass]

	def _createGroups(self, mediaEx, biomass, atpId):
		""" Create metabolite and reaction groups that we'll need (e.g., to set up the S matrix) """

		self.metGroupNew("real", [x for x in self._metIds if x[:len("pseudo")] != "pseudo"])
		self.metGroupNew("pseudo_(f_atp-f_i)", [x for x in self._metIds if x[:len("pseudo_(f_atp-")] == "pseudo_(f_atp-"])
		self.metGroupNew("pseudo_(g_bio-f_i)", [x for x in self._metIds if x[:len("pseudo_(g_bio-")] == "pseudo_(g_bio-"])
		self.metGroupNew("mediaEx", [x["met"] for x in mediaEx])
		self.metGroupNew("biomass", [x["id"] for x in biomass])
		# self.metGroupNew("biomass_not_atp", [x["id"] for x in biomass if x["id"] != atpId])
		# self.metGroupNew("atp", [atpId])

		self.rxnGroupNew("real", [x for x in self._rxnIds if x[:len("rnx_")] == "rxn_"])
		self.rxnGroupNew("mediaEx", [x for x in self._rxnIds if x[:len("mediaEx_")] == "mediaEx_"])
		self.rxnGroupNew("f", [x for x in self._rxnIds if x[:len("r_")] == "f_"])
		self.rxnGroupNew("x", [x for x in self._rxnIds if x[:len("x_")] == "x_"])
		self.rxnGroupNew("(f_atp-f_i)", [x for x in self._rxnIds if x[:len("(f_atp-")] == "(f_atp-"])
		self.rxnGroupNew("g_bio", ["g_bio"])
		self.rxnGroupNew("(g_bio-f_i)", [x for x in self._rxnIds if x[:len("(g_bio-f_")] == "(g_bio-f_"])
		self.rxnGroupNew("f_atp", ["f_" + atpId])
		self.rxnGroupNew("f_not_atp", [x for x in self.rxnGroup("f").ids() if x != "f_" + atpId])
		self.rxnGroupNew("lowerMutable",
						self.rxnGroup("real").ids() + 
						self.rxnGroup("mediaEx").ids() +
						self.rxnGroup("(g_bio-f_i)").ids() +
						self.rxnGroup("x").ids())
		self.rxnGroupNew("lowerImmutable",
						self.rxnGroup("f").ids() +
						self.rxnGroup("(f_atp-f_i)").ids() +
						self.rxnGroup("g_bio").ids())
		self.rxnGroupNew("upperMutable",
						self.rxnGroup("real").ids() + 
						self.rxnGroup("mediaEx").ids() +
						self.rxnGroup("f").ids() +
						self.rxnGroup("(f_atp-f_i)").ids() +
						self.rxnGroup("g_bio").ids())
		self.rxnGroupNew("upperImmutable",
						self.rxnGroup("x").ids() +
						self.rxnGroup("(g_bio-f_i)").ids())

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

		# Populate r reactions
		self._S[self.metGroup("biomass").idxs(), self.rxnGroup("f").idxs()] = np.array([b["coeff"] for b in biomass])

		# Populate x reactions
		self._S[self.metGroup("biomass").idxs(), self.rxnGroup("x").idxs()] = -1.

		# Populate (atp-f_i) reactions
		self._S[self.metGroup("pseudo_(f_atp-f_i)").idxs(), self.rxnGroup("f_atp").idxs()] = 1.
		self._S[self.metGroup("pseudo_(f_atp-f_i)").idxs(), self.rxnGroup("f_not_atp").idxs()] = -1.
		self._S[self.metGroup("pseudo_(f_atp-f_i)").idxs(), self.rxnGroup("(f_atp-f_i)").idxs()] = -1.

		# Don't populate g_bio reaction (it should be all zeros)

		# Populate (g_bio-f_i) reactions
		self._S[self.metGroup("pseudo_(g_bio-f_i)").idxs(), self.rxnGroup("g_bio").idxs()] = 1.
		self._S[self.metGroup("pseudo_(g_bio-f_i)").idxs(), self.rxnGroup("f").idxs()] = -1.
		self._S[self.metGroup("pseudo_(g_bio-f_i)").idxs(), self.rxnGroup("(g_bio-f_i)").idxs()] = -1.

		# S matrix has changed, so any LP solution needs to be (re)calculated
		self._recalculateSolution = True

	def _populateC(self):
		self._c[self.rxnGroup("f_atp").idxs()] = self._alpha
		self._c[self.rxnGroup("g_bio").idxs()] = self._beta
		self._c[self.rxnGroup("(f_atp-f_i)").idxs()] = self._gamma

		self._recalculateSolution = True

	def _populateBounds(self):
		self.v_lowerIs(self.rxnGroup("f").idxs(), 0, True)
		self.v_lowerIs(self.rxnGroup("(f_atp-f_i)").idxs(), 0, True)
		self.v_lowerIs(self.rxnGroup("g_bio").idxs(), 0, True)
		self.v_upperIs(self.rxnGroup("(g_bio-f_i)").idxs(), 0, True)
		self.v_upperIs(self.rxnGroup("x").idxs(), 0, True)

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

		if not overrideImmutable and not np.all(self._v_upper[self.rxnGroup("upperImmutable").idxs()] == 0):
			#import ipdb; ipdb.set_trace()
			raise Exception, "Attempting to change immutable bound."

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

			v_lower = np.maximum(self._v_lower / self._scaleFactor, -10000)
			v_upper = np.minimum(self._v_upper / self._scaleFactor,  10000)

			self._S[self.metGroup("biomass").idxs(), self.rxnGroup("f").idxs()] /= self._scaleFactor

			self._v, tmp = wholecell.utils.linear_programming.linearProgramming(
			"maximize", self._c, self._S, self._b,
			v_lower, v_upper, self._metConstTypes, self._rxnVarTypes,
			None	# Will use glpk
			)

			self._S[self.metGroup("biomass").idxs(), self.rxnGroup("f").idxs()] *= self._scaleFactor

			self._recalculateSolution = False

		return self._v

	def metaboliteProduction(self, metIdxs, solution):
		# NOTE: Users of this function have to be conscious of units and perform
		# any necessary unit conversions themselves.
		P = np.zeros_like(self._S)
		P[:, self.rxnGroup("mediaEx").idxs()] = self._S[:, self.rxnGroup("mediaEx").idxs()]
		P[:, self.rxnGroup("f").idxs()] = self._S[:, self.rxnGroup("f").idxs()]
		P[:, self.rxnGroup("x").idxs()] = self._S[:, self.rxnGroup("x").idxs()]

		# Sign conventions are fun
		return -np.dot(P, solution)[metIdxs]
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