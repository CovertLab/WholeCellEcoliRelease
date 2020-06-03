from __future__ import absolute_import, division, print_function

import numpy as np

from wholecell.utils.modular_fba import FluxBalanceAnalysis

KCAT_MAX = 1.4e6

TOY_MODEL_REACTION_STOICH = {
	"R1": {"A":-1, "ATP":-1, "B":1},
	"R2a": {"B":-1, "ATP":2, "NADH":2, "C":1},
	"R2b": {"C":-1, "ATP":-2, "NADH":-2, "B":1},
	"R3": {"B":-1, "F":1},
	"R4": {"C":-1, "G":1},
	"R5": {"G":-1, "C":0.8, "NADH":2},
	"R6": {"C":-1, "ATP":2, "D":3},
	"R7": {"C":-1, "NADH":-4, "E":3},
	"R8a": {"G":-1, "ATP":-1, "NADH":-2, "H":1},
	"R8b": {"G":1, "ATP":1, "NADH":2, "H":-1},
	"Rres": {"NADH":-1, "O2":-1, "ATP":1},
}

BIOMASS_REACTION_STOICH = {
	"v_biomass": {"C":1, "F":1, "H":1, "ATP":10}
}

TRANSPORT_LIMITS = {
	"A": 21.,
	"F": 5.0,
	"D": -12.0,
	"E": -12.0,
	"H": 5.0,
	"O2": 15.0,
}

REACTION_ENZYMES = {
	"R1":"E1",
	"R2a":"E2a",
	"R2b":"E2b",
	"R3":"E3",
	"R4":"E4",
	"R5":"E5",
	"R6":"E6",
	"R7":"E7",
	"R8a":"E8a",
	"R8b":"E8b",
	"Rres":"Eres",
}

ENZYME_CONCENTRATIONS = {
	"E1":10.0,
	"E2a":10.0,
	"E2b":10.0,
	"E3":10.0,
	"E4":10.0,
	"E5":10.0,
	"E6":10.0,
	"E7":10.0,
	"E8a":10.0,
	"E8b":10.0,
	"Eres":10.0,
}


def testModel(
		toyModelReactionStoich=None,
		biomassReactionStoich=None,
		transportLimits=None,
		reactionEnzymes=None,
		enzymeConcentrations=None):
	if toyModelReactionStoich is None:
		toyModelReactionStoich = TOY_MODEL_REACTION_STOICH
	if biomassReactionStoich is None:
		biomassReactionStoich = BIOMASS_REACTION_STOICH
	if transportLimits is None:
		transportLimits = TRANSPORT_LIMITS
	if reactionEnzymes is None:
		reactionEnzymes = REACTION_ENZYMES
	if enzymeConcentrations is None:
		enzymeConcentrations = ENZYME_CONCENTRATIONS

	fba = FluxBalanceAnalysis(
		reactionStoich=toyModelReactionStoich,
		externalExchangedMolecules=transportLimits.keys(),
		objective=biomassReactionStoich["v_biomass"],
		objectiveType="standard",
		solver="glpk")
	exchangeMolecules = fba.getExternalMoleculeIDs()
	fba.setExternalMoleculeLevels([transportLimits[molID] for molID in exchangeMolecules])

	for reactionID in toyModelReactionStoich:
		if reactionID in reactionEnzymes:
			enzymeID = reactionEnzymes[reactionID]
			if enzymeID in enzymeConcentrations:
				fba.setMaxReactionFlux(reactionID, KCAT_MAX * enzymeConcentrations[enzymeID])

	return fba.getBiomassReactionFlux()[0]

unconstrainedFlux = testModel()

biomassRatesDict = {}
for enzymeID in ENZYME_CONCENTRATIONS:
	testConcentrations = ENZYME_CONCENTRATIONS.copy()
	testConcentrations[enzymeID] = 0.
	rate = testModel(enzymeConcentrations=testConcentrations)
	biomassRatesDict[enzymeID] = rate
	for enzyme2ID in ENZYME_CONCENTRATIONS:
		if enzymeID == enzyme2ID or '('+enzyme2ID + ',' + enzymeID+')' in biomassRatesDict:
			continue
		testConcentrations[enzyme2ID] = 0.
		rate = testModel(enzymeConcentrations=testConcentrations)
		biomassRatesDict['('+enzymeID + ',' + enzyme2ID+')'] = rate


print("Wildtype flux is {}.".format(unconstrainedFlux))
for enzyme in sorted(biomassRatesDict):
	rate = biomassRatesDict[enzyme]
	if np.abs(rate - unconstrainedFlux) > 1e-4:
		print("{} affects growth, when knocked out the biomass reaction flux is {}.".format(enzyme, rate))
