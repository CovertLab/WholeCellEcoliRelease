import numpy as np

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

transportReactionStoich = {
	"Tc1": {"A":1},
	"Tf": {"F":1},
	"Td": {"D":-1},
	"Te": {"E":-1},
	"Th": {"H":1},
	"To2": {"O2":1},
}

biomassReactionStoich = {
	"v_biomass": {"C":1, "F":1, "H":1, "ATP":10}
}

transportLimits = {
	"Tc1": 10.5,
	"Tf": 5.0,
	"Td": 12.0,
	"Te": 12.0,
	"Th": 5.0,
	"To2": 15.0,
}

metaboliteToTransportReaction = {metaboliteID:rxnID for rxnID, rxnStoich in transportReactionStoich.iteritems() for metaboliteID in rxnStoich}

maxReactionFlux = {metaboliteID:(coeff*transportLimits[rxnID]) for rxnID, reactionStoich in transportReactionStoich.iteritems() for metaboliteID, coeff in reactionStoich.iteritems()}

fba = FluxBalanceAnalysis(
	reactionStoich=toyModelReactionStoich,
	externalExchangedMolecules=[molID for rxnStoich in transportReactionStoich.values() for molID in rxnStoich.keys()],
	objective=biomassReactionStoich.values()[0],
	objectiveType="standard",
	)

exchangeMolecules = fba.externalMoleculeIDs()

exchangeLimits = [maxReactionFlux[x] for x in exchangeMolecules]