from __future__ import absolute_import, division, print_function

from wholecell.utils.modular_fba import FluxBalanceAnalysis

KCAT_MAX = 1.4e6

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

toyModelReactionStoichWithBiomass = {
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
	"v_biomass": {"C":-1, "F":-1, "H":-1, "ATP":-10}
}

biomassReactionStoich = {
	"v_biomass": {"C":1, "F":1, "H":1, "ATP":10}
}

transportLimits = {
	"A": 21.,
	"F": 5.0,
	"D": -12.0,
	"E": -12.0,
	"H": 5.0,
	"O2": 15.0,
}

reactionEnzymes = {
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

enzymeConcentrations = {
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

enzymeKcats = {
	"E1":1.0,
	"E2a":.5,
	"E2b":.5,
	"E3":.01,
	"E4":.01,
	"E5":.2,
	"E6":.3,
	"E7":.4,
	"E8a":.4,
	"E8b":.4,
	"Eres":.2,
}

def checkErrors(targetFluxes,
		fixedReactionNames=None,
		reactionStoichiometry=None,
		transportLimits_=None):
	if fixedReactionNames is None:
		fixedReactionNames = ["v_biomass"]
	if reactionStoichiometry is None:
		reactionStoichiometry = toyModelReactionStoichWithBiomass
	if transportLimits_ is None:
		transportLimits_ = transportLimits

	fba_moma = FluxBalanceAnalysis(
		reactionStoich=reactionStoichiometry,
		externalExchangedMolecules=transportLimits_.keys(),
		objective=targetFluxes,
		objectiveType="moma",
		objectiveParameters={"fixedReactionNames":fixedReactionNames},
		solver="glpk")
	exchangeMolecules = fba_moma.getExternalMoleculeIDs()
	fba_moma.setExternalMoleculeLevels([transportLimits_[molID] for molID in exchangeMolecules])
	return fba_moma.errorFluxes(), fba_moma.errorAdjustedReactionFluxes()

fba = FluxBalanceAnalysis(
	reactionStoich=toyModelReactionStoich,
	externalExchangedMolecules=transportLimits.keys(),
	objective=biomassReactionStoich["v_biomass"],
	objectiveType="standard",
	solver="glpk",
)
exchangeMolecules = fba.getExternalMoleculeIDs()
fba.setExternalMoleculeLevels([transportLimits[molID] for molID in exchangeMolecules])
wildtypeBiomassFlux = fba.getBiomassReactionFlux()

# Adjust kcats
targetFluxes = {}
for reactionID, enzymeID in reactionEnzymes.iteritems():
	targetFluxes[reactionID] = enzymeConcentrations[enzymeID] * enzymeKcats[enzymeID]
targetFluxes["v_biomass"] = wildtypeBiomassFlux

errors, rates = checkErrors(targetFluxes)

errors_dict = dict(zip(enzymeKcats, errors))

kcat_adjustments = {enzymeID: error / enzymeConcentrations[enzymeID] for enzymeID, error in errors_dict.iteritems()}

for enzymeID, error in kcat_adjustments.iteritems():
	print("{} kcat error is {}.".format(enzymeID, error))
