from __future__ import absolute_import, division, print_function

from wholecell.utils.modular_fba import FluxBalanceAnalysis

toyModelReactionStoich = {
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

biomassReactionStoich = {
	"v_biomass": {"C":1, "F":1, "H":1, "ATP":10}
}

transportLimits = {
	"A": 21.0,
	"F": 5.0,
	"D": -12.0,
	"E": -12.0,
	"H": 5.0,
	"O2": 15.0,
}


fba = FluxBalanceAnalysis(
	reactionStoich=toyModelReactionStoich,
	externalExchangedMolecules=list(transportLimits.keys()),
	objective=biomassReactionStoich["v_biomass"],
	objectiveType="standard",
	solver="glpk",
	)

exchangeMolecules = fba.getExternalMoleculeIDs()

fba.setExternalMoleculeLevels([transportLimits[molID] for molID in exchangeMolecules])

biomassReactionFlux = fba.getBiomassReactionFlux()[0]

print(biomassReactionFlux)
