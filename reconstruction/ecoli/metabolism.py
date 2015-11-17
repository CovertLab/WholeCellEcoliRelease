
from __future__ import division

import collections

import numpy as np

from wholecell.utils import units
from hard_coded_data import *

class Metabolism(object):
	def __init__(self, kb):
		# Build the matrices/vectors for metabolism (FBA)

		# These may be modified/extended later, but should provide the basic
		# data structures

		exchangeReactions = [x["id"] for x in kb._reactions if x["id"].startswith("FEIST_EX") or x["id"].startswith("FEIST_DM_")]
		exchangeReactions += ["SELNP_MEDIA_EXCHANGE_HACKED"]

		disabledReactions = ("FEIST_CAT_0", "FEIST_SPODM_0", "FEIST_SPODMpp",
			"FEIST_FHL_0_0", "FEIST_FHL_1_0", "FEIST_ATPM")
		# NOTE: the first five were disabled by Feist, the last is NGAM which model explictly elsewhere

		reactionStoich = {}
		externalExchangeMolecules = []
		reversibleReactions = []
		reactionEnzymes = {}
		reactionRates = {}

		unconstrainedExchangeMolecules = ("CA2[e]", "CL[e]", "CO2[e]", "COBALT2[e]",
			"CU2[e]", "FE2[e]", "FE3[e]", "H[e]", "H2O[e]", "K[e]", "MG2[e]",
			"MN2[e]", "MOBD[e]", "NA1[e]", "NH4[e]", "PI[e]", "SO4[e]", "TUNGS[e]",
			"ZN2[e]", "SELNP[e]",)

		exchangeUnits = units.mmol / units.g / units.h

		constrainedExchangeMolecules = {
			"CBL1[e]": exchangeUnits * 0.01, # TODO: try removing this constraint
			"GLC-D[e]": exchangeUnits * 8,
			"O2[e]": exchangeUnits * 18.5
			}

		catalysisUnits = 1 / units.s

		validEnzymeIDs = set(kb.bulkMolecules["moleculeId"])
		validEnzymeCompartments = collections.defaultdict(set)

		for enzymeID in validEnzymeIDs:
			enzyme = enzymeID[:enzymeID.index("[")]
			location = enzymeID[enzymeID.index("[")+1:enzymeID.index("[")+2]

			validEnzymeCompartments[enzyme].add(location)

		for reaction in kb._reactions:
			reactionID = reaction["id"]
			stoich = reaction["stoichiometry"]

			if reactionID in disabledReactions:
				continue

			elif reactionID in exchangeReactions:
				if len(stoich) != 1:
					raise Exception("Invalid exchange reaction")

				externalExchangeMolecules.append("{}[{}]".format(
					stoich[0]["molecule"], stoich[0]["location"]
					))

			else:
				if len(stoich) <= 1:
					raise Exception("Invalid biochemical reaction")

				reducedStoich = {
					"{}[{}]".format(
						entry["molecule"], entry["location"]
						): entry["coeff"]
					for entry in stoich
					}

				reactionStoich[reactionID] = reducedStoich

				# Assign reversibilty
				if reaction["dir"] == 0:
					reversibleReactions.append(reactionID)

				# Assign k_cat, if known
				kcat = reaction["kcat"]

				if kcat is not None:
					reactionRates[reactionID] = kcat * catalysisUnits

				# Assign enzyme, if any
				if reactionID in REACTION_ENZYME_ASSOCIATIONS.viewkeys():
					enzymes = REACTION_ENZYME_ASSOCIATIONS[reactionID]

				else:
					enzymes = reaction['catBy']

				reactantLocations = {reactant["location"]
					for reactant in reaction["stoichiometry"]}

				if enzymes is not None and len(enzymes) > 0:
					if len(enzymes) > 1:
						raise Exception("Reaction {} has multiple associated enzymes: {}".format(
							reactionID, enzymes))

					(enzyme,) = enzymes

					validLocations = validEnzymeCompartments[enzyme]

					if len(validLocations) == 0:
						raise Exception("Reaction {} uses enzyme {} but this enzyme does not exist.".format(
							reactionID,
							enzyme
							))

					if len(validLocations) == 1:
						(location,) = validLocations

					elif len(reactantLocations) == 1:
						(location,) = reactantLocations

					elif reactantLocations == {"p", "e"}: # if reaction is periplasm <-> extracellular
						(location,) = {"o"} # assume enzyme is in outer membrane

					elif reactantLocations == {"c", "p"}: # if reaction is cytoplasm <-> periplasm
						(location,) = {"i"} # assume enzyme is in inner membrane

					else:
						raise Exception("Reaction {} has multiple associated locations: {}".format(
							reactionID,
							reactantLocations
							))

					assert location in validLocations

					enzymeID = "{}[{}]".format(enzyme, location)

					reactionEnzymes[reactionID] = enzymeID

		mws = kb.getter.getMass(externalExchangeMolecules)

		exchangeMasses = {moleculeID:mws[index]
			for index, moleculeID in enumerate(externalExchangeMolecules)}

		# Filter out reaction-enzyme associations that lack rates

		reactionEnzymes = {
			reactionID:enzymeID
			for reactionID, enzymeID in reactionEnzymes.viewitems()
			if reactionRates.has_key(reactionID)
			}

		self.reactionStoich = reactionStoich
		self.externalExchangeMolecules = externalExchangeMolecules
		self._exchangeMasses = exchangeMasses
		self.reversibleReactions = reversibleReactions
		self.reactionEnzymes = reactionEnzymes
		self._reactionRates = reactionRates
		self._unconstrainedExchangeMolecules = unconstrainedExchangeMolecules
		self._constrainedExchangeMolecules = constrainedExchangeMolecules


	def reactionRates(self, timeStep):
		return {
			reactionID:(reactionRate * timeStep).asNumber()
			for reactionID, reactionRate in self._reactionRates.viewitems()
			}


	def exchangeMasses(self, targetUnits):
		return {
			moleculeID:mass.asNumber(targetUnits)
			for moleculeID, mass in self._exchangeMasses.viewitems()
			}


	def exchangeConstraints(self, exchangeIDs, coefficient, targetUnits):
		externalMoleculeLevels = np.zeros(len(exchangeIDs), np.float64)

		for index, moleculeID in enumerate(exchangeIDs):
			if moleculeID in self._unconstrainedExchangeMolecules:
				externalMoleculeLevels[index] = np.inf

			elif moleculeID in self._constrainedExchangeMolecules.viewkeys():
				externalMoleculeLevels[index] = (
					self._constrainedExchangeMolecules[moleculeID] * coefficient
					).asNumber(targetUnits)

		return externalMoleculeLevels
