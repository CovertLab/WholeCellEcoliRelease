#!/usr/bin/env python

import numpy as np

from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli

from wholecell.utils import units

#equation string interpreter
from Equation import Expression




# class EnzymeKinetics(wholecell.processes.process.Process):
# 	""" EnzymeKinetics """

# 	_name = "EnzymeKinetics"

# 	# Constructor
# 	def __init__(self):
# 		# Views

# 		self.metaboliteConcentrations = None
# 		self.enzymeConcentrations = None

# 		super(EnzymeKinetics, self).__init__()


# 	# Construct object graph
# 	def initialize(self, sim, kb):
# 		super(EnzymeKinetics, self).initialize(sim, kb)

# 		# Load reaction rate law information

# 		reactions = kb.process.metabolism.reactionRateInfo

# 		# Views

# 		# SET VIEWS!

# 		self.metaboliteConcentrations = self.bulkMoleculeView(kb.moleculeGroups.s30_fullComplex[0])
# 		self.ribosome50S = self.bulkMoleculeView(kb.moleculeGroups.s50_fullComplex[0])

# 		self.mRnas = self.bulkMoleculesView(mrnaIds)










# Returns the michaelis-Menton predicted rate of a reaction.
# Inputs: 	enzyme_conc (concentration of the enzyme)
# 			substrate_conc (concentration of the reaction substrate/metabolite)
# 			k_cat (max catalytic rate of the enzyme)
# 			k_M (Michaelis-Menton constant for the enzyme (substrate concentration at which half of max rate is reached))
# Returns: 	predicted rate of the reaction
def michaelisMenton(enzyme_conc, substrate_conc, k_cat, k_M):
	return (k_cat*enzyme_conc)*((substrate_conc)/(k_M + substrate_conc))


# Returns the theoretical maximum catalytic rate of an enzymatic reaction, for use in jFBA as an upper bound.
# Inputs: 	enzyme_conc (concentration of the enzyme)
# 			k_cat (max catalytic rate of the enzyme)
# Returns: 	Maximum theoretical rate of the reaction
def maxReactionRate(enzyme_conc, k_cat):
	return k_cat*enzyme_conc


# Returns the approximated rate of a reaction with 1 or more substrates and 0 or more inhibitors.
# Inputs: 	enzyme_conc (concentration of the enzyme)
# 			k_cat (max catalytic rate of the enzyme)
# 			substrate_conc_array (array of concentrations of each reaction substrate/metabolite)
# 			k_M_array (array of Michaelis-Menton constants for the enzyme with each substrate)
# 			inhibitor_conc_array (array of concentrations of each reaction inhibitor)
# 			k_I_array (array of Michaelis-Menton constants for the enzyme with each inhibitor)
# Returns: 	predicted rate of the reaction.
# Notes: 	Assumes that substrate_conc_array and inhibitor_conc_array are 1-to-1 with k_M_array
# 			and k_I_array respectively, and are in the same order.
# 			If no inhibitor is to be used, input an empty array for inhibitor_conc_array and k_I_array.
#			If no inputs given for inhibitor_conc_array and k_I_array, function defaults to this behavior.
# 			Must have at least one substrate, throws an error if no substrates in substrate_conc_array
def enzymeRateApproximate(enzyme_conc, k_cat, substrate_conc_array, k_M_array, inhibitor_conc_array = [], k_I_array = []):
	rate = enzyme_conc*k_cat
	
	#  Require that at least one substrate be used
	assert (len(substrate_conc_array) > 0)

	# Check that the same number of substrates and Michaelis-Menton constants are given
	assert (len(substrate_conc_array) == len(k_M_array))

	# Check that the same number of inhibitors and inhibitor Michaelis-Menton constants are given
	assert (len(inhibitor_conc_array) == len(k_I_array))

	# Adjust rate for all substrates
	for i, k_M in enumerate(k_M_array):
		rate *= ((substrate_conc_array[i])/(k_M + substrate_conc_array[i]))

	# Adjust rate for all inhibitors
	for i, k_I in enumerate(k_I_array):
		rate *= ((inhibitor_conc_array[i])/(k_I + inhibitor_conc_array[i]))
	
	return rate	


# Returns the approximated rate of a reaction using a rate equation passed in as a string.
# Inputs: 	enzyme_conc (concentration of the enzyme)
# 			k_cat (max catalytic rate of the enzyme)
# 			eq_string (a string corresponding to the rate equation of the
# 						reaction - see notes for formatting)
# 			parameter_definition_array (an array of strings of the variables appearing in eq_string)
# 			parameters_array (an array of the parameters to be plugged in to the rate 
# 							equation, in the same order in which they appear in parameter_definition_array)
# Returns: 	predicted rate of the reaction.
# Notes: 	This function should only be used if none of the others are able to define the rate equation!
# 				In particular, enzymeRateApproximate() can flexibly define many rate laws, try that one first.
# 			Formatting of eq_string is as a valid python mathematical expression, for example:
# 				"k_cat*E*((substrate_conc)/(substrate_conc + k_M))"
# 			Formatting of parameter_definition_array must match the variables in eq_string, eg:
# 				"["E","k_cat","substrate_conc","k_M"]"
# 			Then parameters_array must match the order of parameter_definition_array, eg:
# 				"["2.0",".0343","3.9",".7"]" --> E is 2.0, k_cat is .0343, etc
# 			Function checks that there are equal numbers of parameters and parameter definitions, but
# 				not that they are in the right order.
# 			Function does NOT check for validity of the equation passed in in eq_string.
def enzymeRateCustom(eq_string, parameter_definition_array, parameters_array):
	# Check that there are equal numbers of parameter values and parameter definitions
	assert (len(parameters_array) == len(parameter_definition_array))

	# Set up the custom function
	customRateLaw = Expression(eq_string, parameter_definition_array)

	# Evaluate the custom function
	return customRateLaw(*parameters_array)


def enzymeRate(reactionInfo, enzymeConcArray, substrateConcArray, inhibitorConcArray = []):
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
	"""
	# reactionInfo = kb.process.metabolism.reactionRateInfo[reactionID]

	# Standard or custom reaction rate law?
	if(reactionInfo["rateEquationType"] == "standard"):
		# Standard rate law
		rate = enzymeRateApproximate(enzymeConcArray[0], reactionInfo["kcat"], substrateConcArray, reactionInfo["kM"], inhibitorConcArray, reactionInfo["kI"])
	
	elif(reactionInfo["rateEquationType"] == "custom"):
		# Custom rate law
		equationString = reactionInfo["customRateEquation"]
		parameterDefinitionArray = reactionInfo["customParameters"]
		parametersArray = reactionInfo["customParameterConstantValues"] + enzymeConcArray + substrateConcArray

		rate = enzymeRateCustom(equationString, parameterDefinitionArray, parametersArray)

	return rate