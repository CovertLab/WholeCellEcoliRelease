#!/usr/bin/env python

import numpy as np

from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli

from wholecell.utils import units

import theano.tensor as T

# Equation string interpreter
from Equation import Expression




# class EnzymeKinetics(object):
# 	"""
# 	EnzymeKinetics

# 	Stores a compiled theano function determining any reaction kinetics known
# 	"""



# 	def __init__(self, kb):
# 		self.reactions = [] # sorted list of reactions with known kinetics
# 		self.rateFunctions = [] # sorted list of rate functions

# 		# Load rate functions from enzymeKinetics.tsv flat file
# 		self.reactionRateInfo = kb.process.metabolism.reactionRateInfo
# 		self.enzymesWithKineticInfo = kb.process.metabolism.enzymesWithKineticInfo["enzymes"]

# 		# Load info on all reactions in the model
# 		self.allReactions = kb.process.metabolism.reactionStoich


# 		# find reaction rate limits

# 		# This list will hold the limits - default limit to infinity
# 		self.reactionRates = np.ones(len(self.fba.reactionIDs()))*np.inf
# 		self.perEnzymeRates = np.ones(len(self.fba.reactionIDs()))*np.inf
# 		self.enzymeConc = np.zeros(len(self.fba.reactionIDs()))

# 		for index, reactionID in enumerate(self.fba.reactionIDs()):
# 			rateInfo = {}
# 			try:
# 				rateInfo = self.reactionRateInfo[reactionID]
# 			except:
# 				continue

# 			substrateIDs = rateInfo["substrateIDs"]
# 			substrateIXs = [self.metaboliteIndexDict[x] for x in substrateIDs]
# 			substrateConcArray = [metaboliteConcentrations[i] for i in substrateIXs]

# 			enzymeIDs = rateInfo["enzymeIDs"]
# 			enzymeIXs = [self.enzymeIndexDict[x] for x in enzymeIDs]
# 			enzymeConcArray = [enzymeConcentrations[i] for i in enzymeIXs]

# 			rate = enzymeRate(rateInfo, enzymeConcArray, substrateConcArray)

# 			self.reactionRates[index] = rate

# 			# Assumes the least concentrated enzyme limit rate, if multiple
# 			# (Almost always this will be the single enzyme in the rate law)
# 			self.perEnzymeRates[index] = rate / np.amin(enzymeConcArray)

# 			self.enzymeConc[index] = np.amin(enzymeConcArray)



def michaelisMenton(enzyme_conc, substrate_conc, k_cat, k_M):
	""" Returns the michaelis-Menton predicted rate of a reaction.
		
		Inputs: enzyme_conc (concentration of the enzyme)
				substrate_conc (concentration of the reaction substrate/metabolite)
				k_cat (max catalytic rate of the enzyme)
				k_M (Michaelis-Menton constant for the enzyme (substrate concentration at which half of max rate is reached))
		Returns: 	predicted rate of the reaction
	"""
	
	return (k_cat*enzyme_conc)*((substrate_conc)/(k_M + substrate_conc))


def maxReactionRate(enzyme_conc, k_cat):
	""" Returns the theoretical maximum catalytic rate of an enzymatic reaction, for use in jFBA as an upper bound.
		
		Inputs:  enzyme_conc (concentration of the enzyme)
				 k_cat (max catalytic rate of the enzyme)
		Returns: Maximum theoretical rate of the reaction
	"""
	
	return k_cat*enzyme_conc


def enzymeRateApproximate(enzyme_conc, k_cat, substrate_conc_array, k_M_array, k_I_array = []):

	""" Returns the approximated rate of a reaction with 1 or more substrates and 0 or more inhibitors.
		
		Inputs: 	enzyme_conc (concentration of the enzyme)
			k_cat (max catalytic rate of the enzyme)
			substrate_conc_array (array of concentrations of each reaction substrate/metabolite/inhibitor)
			k_M_array (array of Michaelis-Menton constants for the enzyme with each substrate)
			k_I_array (array of Michaelis-Menton constants for the enzyme with each inhibitor)
		Returns: 	predicted rate of the reaction.
		Notes: 	Assumes that substrate_conc_array is in the same order as k_M_array followed by k_I_array.
			In other words, if two k_M's are given and one k_I, then the substrate_conc_array must be in
			the order corresponding to: [k_M_1, k_M_2, k_I_1]. 
			If no inhibitor is to be used, input an empty array for k_I_array.
			Must have at least one substrate, throws an error if no substrates in substrate_conc_array.
			Inhibitors are noncompetitive in this function.
	"""
	
	rate = enzyme_conc*k_cat

	if (len(substrate_conc_array) != len(k_M_array)):
		import ipdb; ipdb.set_trace()
	
	#  Require that at least one substrate be used
	assert (len(substrate_conc_array) > 0)

	# Check that the same number of substrates and Michaelis-Menton constants are given
	assert (len(substrate_conc_array) == len(k_M_array) + len(k_I_array))

	n = 0
	# Adjust rate for all substrates
	for k_M in k_M_array:
		rate *= ((substrate_conc_array[n])/(k_M + substrate_conc_array[n]))
		n += 1

	# Adjust rate for any/all inhibitors
	for k_I in k_I_array:
		rate *= ((1)/(1 + (substrate_conc_array[n]/k_I)))
		n += 1

	return rate	

def enzymeRateCustom(eq_string, parameter_definition_array, parameters_array):
	""" Returns the approximated rate of a reaction using a rate equation passed in as a string.
		
		Inputs: 	enzyme_conc (concentration of the enzyme)
				k_cat (max catalytic rate of the enzyme)
				eq_string (a string corresponding to the rate equation of the
							reaction - see notes for formatting)
				parameter_definition_array (an array of strings of the variables appearing in eq_string)
				parameters_array (an array of the parameters to be plugged in to the rate 
								equation, in the same order in which they appear in parameter_definition_array)
		Returns: 	predicted rate of the reaction.
		Notes: 	This function should only be used if none of the others are able to define the rate equation!
					In particular, enzymeRateApproximate() can flexibly define many rate laws, try that one first.
				Formatting of eq_string is as a valid python mathematical expression, for example:
					"k_cat*E*((substrate_conc)/(substrate_conc + k_M))"
				Formatting of parameter_definition_array must match the variables in eq_string, eg:
					"["E","k_cat","substrate_conc","k_M"]"
				Then parameters_array must match the order of parameter_definition_array, eg:
					"["2.0",".0343","3.9",".7"]" --> E is 2.0, k_cat is .0343, etc
				Function checks that there are equal numbers of parameters and parameter definitions, but
					not that they are in the right order.
			Function does NOT check for validity of the equation passed in in eq_string.
	"""

	import ipdb; ipdb.set_trace()

	# Check that there are equal numbers of parameter values and parameter definitions
	assert (len(parameters_array) == len(parameter_definition_array))

	# Set up the custom function
	customRateLaw = Expression(eq_string, parameter_definition_array)

	# Evaluate the custom function
	return customRateLaw(*parameters_array)


def enzymeRate(reactionInfo, enzymeConcArray, substrateConcArray):
	""" Returns the approximated rate of a reaction.
		
		Inputs: 	reactionInfo (Dict of info for reaction whose rate limit is
						to be calculated)
					enzymeConcArray (array of concentrations of enzymes)
						NOTE: if using custom rate law: these must be in the
						same order as they appear in the "customParameters"
						array of the reactionInfo dictionary
					substrateConcArray (array of concentrations of substrates)
						NOTE: if using custom rate law: these must be in the
						same order as they appear in the "customParameters"
						array of the reactionInfo dictionary
					inhibitorConcArray (array of concentrations of inhibitors)
						NOTE: Only us for standard rate law. Must be in the
						same order as their corresponding kIs in the in the
						"kI" array of the reactionInfo dictionary
		Returns: predicted rate of the reaction.
		Notes:	
				This function has very different behavior for custom and for
					standard rate laws (determined by the "rateEquationType"
					field of reactionInfo)
				For custom rate laws: the order of
					reactionInfo["customParameters"] is assumed to be enzyme
					concentrations followed by substrate concentrations (in the
					same order passed in as substrateConcArray.
				For standard rate laws: the reaction is assumed to have only
					one enzyme, and only the first enzyme in enzymeConcArray
					is considered.
				Currently uses ONLY the higest kcat value, if multiple are available.
	"""

	# Standard or custom reaction rate law?
	if(reactionInfo["rateEquationType"] == "standard"):
		# Standard rate law

		# Check if K_M or K_I given
		if(len(reactionInfo["kM"]) + len(reactionInfo["kI"])>1):		
			# Use michaelis-menton kinetics
			rate = enzymeRateApproximate(enzymeConcArray[0], np.amax(reactionInfo["kcat"]), substrateConcArray, reactionInfo["kM"], reactionInfo["kI"])
		else:
			# Use only the kcat
			rate = maxReactionRate(enzymeConcArray[0], np.amax(reactionInfo["kcat"]))
	elif(reactionInfo["rateEquationType"] == "custom"):
		# Custom rate law
		equationString = reactionInfo["customRateEquation"]
		parameterDefinitionArray = reactionInfo["customParameters"]
		parametersArray = reactionInfo["customParameterConstantValues"] + enzymeConcArray + substrateConcArray

		rate = enzymeRateCustom(equationString, parameterDefinitionArray, parametersArray)

	return rate