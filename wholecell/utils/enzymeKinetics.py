#!/usr/bin/env python

"""
EnzymeKinetics

Compiles a theano function which can be called to determine the rates of all reactions in a metabolic model.

@author: Morgan Paull
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 8/4/2015
"""

import numpy as np

from wholecell.utils import units

import theano.tensor as T
from theano import function

import re


class EnzymeKinetics(object):
	"""
	EnzymeKinetics

	Stores a compiled theano function determining any reaction kinetics known.
	"""

	def __init__(self, enzymesWithKineticInfo, constraintIDs, reactionRateInfo, reactionIDs, metaboliteIDs, kcatOnly=False):

		# Set default reaction rate limit, to which reactions are set absent other information
		self.defaultRate = np.inf

		# Load rate functions from enzymeKinetics.tsv flat file
		self.reactionRateInfo = reactionRateInfo
		self.constraintIDs = constraintIDs
		self.enzymesWithKineticInfo = enzymesWithKineticInfo

		# Make a dictionary mapping a substrate ID to it's index in self.metabolites()
		self.metaboliteIndexDict = {}
		substrate_vars_array = [0]*len(metaboliteIDs)
		for index,name in enumerate(metaboliteIDs):
			self.metaboliteIndexDict[name] = index
			substrate_vars_array[index] = T.dscalar(name + '_concentration')

		# Make a dictionary mapping an enzyme ID to it's index in self.enzymes()
		self.enzymeIndexDict = {}
		enzyme_vars_array = [0]*len(self.enzymesWithKineticInfo)
		for index,name in enumerate(self.enzymesWithKineticInfo):
			self.enzymeIndexDict[name] = index
			enzyme_vars_array[index] = T.dscalar(name + '_concentration')

		# Build a function to determine the rate of reactions known to the model
		noRate = T.dscalar('noRate')
		rateExpressionsArray = [noRate]*len(reactionIDs)
		for index, reactionID in enumerate(reactionIDs):
			rateInfo = {}
			try:
				rateInfo = self.reactionRateInfo[reactionID]
			except:
				continue

			rateExpressionsArray[index] = self.buildRateExpression(rateInfo, enzyme_vars_array, substrate_vars_array, self.metaboliteIndexDict, self.enzymeIndexDict, kcatOnly)

		# Build a function to determine the rate of all possible constraints
		longRateExpressionsArray = [noRate]*len(self.reactionRateInfo)
		for index, constraintID in enumerate(self.constraintIDs):
			rateInfo = self.reactionRateInfo[constraintID]
			longRateExpressionsArray[index] = self.buildRateExpression(rateInfo, enzyme_vars_array, substrate_vars_array, self.metaboliteIndexDict, self.enzymeIndexDict, kcatOnly)



		## Compile a theano function for the enzyme kinetics
		# Inputs: the concentrations of every enzyme in enzymesWIthKineticInfo,
		# 			the concentration of every substrate in metaboliteIDs, and 
		# 			then the default rate for reactions without kinetic info,
		# 			in that order
		# Outputs: an array of kinetic rates of reactions ordered as in
		#			reactionIDs
		self.rateFunction = function(enzyme_vars_array + substrate_vars_array + [noRate], T.stack(rateExpressionsArray), on_unused_input='ignore')


		self.allRatesFunction = function(enzyme_vars_array + substrate_vars_array + [noRate], T.stack(longRateExpressionsArray), on_unused_input='ignore')

	def buildRateExpression(self, rateInfo, enzyme_vars_array, substrate_vars_array, metaboliteIndexDict, enzymeIndexDict, kcatOnly):

		# Find the enzyme variable for this reaction.
		# Only uses the first enzyme if there are more than one.
		enzyme_var = enzyme_vars_array[enzymeIndexDict[rateInfo["enzymeIDs"][0]]]

		# Standard or custom reaction rate law?
		if(rateInfo["rateEquationType"] == "standard"):
			# Standard rate law

			# Check if K_M or K_I given
			if(len(rateInfo["kM"]) + len(rateInfo["kI"])>0) and (kcatOnly==False):		
				# Use michaelis-menton kinetics

				# Build a list of substrates vars for this reaction.
				# Input data must be in same order as [k_M's] then [k_I's]
				specific_substrate_vars_array = []
				for substrate_name in rateInfo["substrateIDs"]:
					specific_substrate_vars_array.append(substrate_vars_array[metaboliteIndexDict[substrate_name]])

				# Find the rate function
				rateExpression = self.enzymeRateApproximate(rateInfo, enzyme_var, specific_substrate_vars_array)
			else:
				# Use only the kcat
				rateExpression = self.maxReactionRate(enzyme_var, np.amax(rateInfo["kcat"]))

		elif(rateInfo["rateEquationType"] == "custom"):
			# Custom rate law
			rateExpression = self.enzymeRateCustom(rateInfo, enzyme_vars_array, substrate_vars_array, metaboliteIndexDict, enzymeIndexDict)

		return rateExpression


	def maxReactionRate(self, enzyme_var, k_cat):
		"""
		Returns the theoretical maximum catalytic rate of an enzymatic reaction, to be compiled into a theano function.

		The k_cat is compiled into the function, the enzyme_var is a thenao tensor
			and left as an input to the function.
		"""

		return k_cat * enzyme_var


	def enzymeRateApproximate(self, rateInfo, enzyme_var, substrate_vars_array):

		""" 
		Returns the approximated rate of a reaction with 1 or more substrates and 0 or more inhibitors.

		Inputs: rateInfo - the object defining the enzyme kinetics of this reaction.
				enzyme_vars_array - theano tensor corresponding to the enzyme used
					in this reaction.
				substrate_vars_array - array of theano tensors corresponding to the
					substrates used in this reaction. Must be in the same order as
					the k_M_array and k_I_array defined in rateInfo, and must be
					in order: k_M substrates followed by k_I substrates.

		Returns: an expression for the kinetics of this reaction, which can be
					turned into a theano function later. Uses Michaelis-Menton
					like kinetics, with each substrate either saturating or
					inhibiting with respect to it's k_M or k_I respoectively.
		"""

		k_cat = np.amax(rateInfo["kcat"])

		rate = k_cat * enzyme_var

		n = 0
		# Adjust rate for all substrates
		for k_M in rateInfo["kM"]:
			rate *= ((substrate_vars_array[n])/(k_M + substrate_vars_array[n]))
			n += 1

		# Adjust rate for any/all inhibitors
		for k_I in rateInfo["kI"]:
			rate *= ((1)/(1 + (substrate_vars_array[n]/k_I)))
			n += 1


		return rate

	def enzymeRateCustom(self, rateInfo, enzyme_vars_array, substrate_vars_array, metaboliteIndexDict, enzymeIndexDict):
		"""
		Given an equation string, returns that expression with theano variables substituted in.
		"""

		# Make a dictionary mapping from given user-defined variable to theano var
		D = {}
		placeholder_dict = rateInfo["customParameterVariables"]
		for placeholder in placeholder_dict:
			ID = placeholder_dict[placeholder]
			if ID in metaboliteIndexDict:
				D[placeholder] = substrate_vars_array[metaboliteIndexDict[ID]]
			if ID in enzymeIndexDict:
				D[placeholder] = enzyme_vars_array[enzymeIndexDict[ID]]


		# If the D dict is not the same size as the placeholder dict, then some
		# enzyme or substrate was not recognized by the model.
		if len(D) != len(placeholder_dict):
			raise NameError("One or more enzymes or substrates is not known to the model: %s" % str(placeholder_dict) )

		# Make a dictionary mapping from symbols for constants to their values.
		constants_dict = {}
		symbol_constants = rateInfo["customParameterConstants"]
		constant_values = rateInfo["customParameterConstantValues"]

		for index, value in enumerate(symbol_constants):
			constants_dict[value] = constant_values[index]

		# Build the rate expression
		rate = rateInfo["customRateEquation"]

		for variable in rateInfo["customParameters"]:
			# if it should be a theano variable
			if variable in D:
				rate = re.sub(r"\b%s\b" % variable, "D['" + variable + "']", rate)

			# If it is a constant
			if variable in constants_dict:
				rate = re.sub(r"\b%s\b" % variable, '(' + str(constants_dict[variable]) + ')', rate)

		return eval(rate)