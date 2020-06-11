# Toy model of jFBA using a simplified metabolic network from Fig. 4 and Table 1 of M. Covert's 2001 paper
# Implemented using modular_FBA

from __future__ import absolute_import, division, print_function

import numpy as np
from wholecell.utils.modular_fba import FluxBalanceAnalysis
import matplotlib.pyplot as plt
from wholecell.analysis.analysis_tools import exportFigure
import os
from six.moves import zip

reactions_file = "reactions.tsv"
transportFluxes_file = "transport_fluxes.tsv"
transportRates_file = "transport_rates.tsv"
FBAObjectives_file = "FBA_objectives.tsv"


# Define file parsing
def parseFile(filename):
	out = []
	F = open(filename, mode = "rU")
	
	for line in F:
		line_list = line.split()
		reaction_name = line_list[0]
		molecule_names = []
		molecule_coeffs = []

		for i in np.arange(1, len(line_list) - 1, 2):
			molecule_names.append(line_list[i])
			molecule_coeffs.append(line_list[i + 1])

		out.append({
			"name": reaction_name,
			"stoichiometry": dict(zip(molecule_names, molecule_coeffs)),
			}
		)
	F.close()
	return out


def getValues(filename):
	out = dict()
	F = open(filename, mode = "rU")

	for line in F:
		line_list = line.split()
		molecule = line_list[0]
		rate = line_list[-1]

		out[molecule] = float(rate)

	return out


def runFBA(internalMoleculeLevels):
	fba.setInternalMoleculeLevels(internalMoleculeLevels)
	deltaMolecules = fba.getOutputMoleculeLevelsChange()
	return deltaMolecules


# Get data from files
current_path = os.path.dirname(os.path.abspath(__file__))
reactions = parseFile("%s/%s" % (current_path, reactions_file))
transportFluxes = parseFile("%s/%s" % (current_path, transportFluxes_file))
transportRates = getValues("%s/%s" % (current_path, transportRates_file))
FBAObjectives = getValues("%s/%s" % (current_path, FBAObjectives_file))


# Describe reaction stoichiometry
# - a dict of strings:dicts (reactionID:reaction stoich). Each value in the dict is a dict of molecule ID to stoichiometry pairs
reactionStoich = dict()
for reaction in reactions:
	reactionName = reaction["name"]
	reactionStoich[reactionName] = dict()

	for molecule in reaction["stoichiometry"]:
		reactionStoich[reactionName][molecule] = reaction["stoichiometry"][molecule]


# Describe external exchanged molecules
# - externalExchangedMolecules, an iterable of strings (moleculeIDs). Every provided ID will be set up with an exchange flux.
externalExchangedMolecules = []
for reaction in transportFluxes:
	for molecule in reaction["stoichiometry"].keys():
		externalExchangedMolecules.append(molecule)

nExternalExchangedMolecules = len(externalExchangedMolecules)


# Describe the objective
# - objective, a dict of strings:floats (moleculeID:objective value). The meaning and usage of the objective will vary depending on the formulation of FBA desired.
objective = FBAObjectives


# Set up FBA solver
fba = FluxBalanceAnalysis(
	reactionStoich,
	externalExchangedMolecules,
	objective,
	objectiveType = "homeostatic",
	)


# Describe reaction flux constraints
reactionIDs = list(reactionStoich.keys())
nReactions = len(reactionIDs)
maxReactionFluxes = np.inf * np.ones(nReactions)

for reactionIndex, reaction in enumerate(reactionIDs):

	# Set upper bound constraint according to transport rates
	if reaction in transportRates:
		maxReactionFluxes[reactionIndex] = transportRates[reaction]

fba.setMinReactionFluxes(reactionIDs, np.zeros(nReactions))
fba.setMaxReactionFluxes(reactionIDs, maxReactionFluxes)


# Describe external molecule levels
moleculeIDs = fba.getInternalMoleculeIDs()
nMolecules = len(moleculeIDs)

externalMoleculeLevels = np.inf * np.ones(nExternalExchangedMolecules)
fba.setExternalMoleculeLevels(externalMoleculeLevels)


# Describe internal molecule levels
moleculeCountsInit = np.zeros(nMolecules)
internalMoleculeLevels = moleculeCountsInit


# Perform simulation
nTimesteps = 15
moleculeCounts = [internalMoleculeLevels]

foods = ["ATP"]
foodUnit = 2

for timestep in np.arange(nTimesteps - 1):
	deltaMolecules = runFBA(internalMoleculeLevels)
	moleculeCountsCurrent = internalMoleculeLevels + np.round(deltaMolecules)
	moleculeCounts.append(moleculeCountsCurrent)

	internalMoleculeLevels = moleculeCountsCurrent

	## Eat
	for food in foods:
		internalMoleculeLevels[moleculeIDs.index(food)] -= foodUnit

moleculeCounts = np.array(moleculeCounts)


# Plot
rows = 4
cols = 3

fig = plt.figure()

for index, molecule in enumerate(moleculeIDs):
	ax = plt.subplot(rows, cols, index + 1)
	ax.plot(np.arange(nTimesteps), moleculeCounts[:, index])

	ax.plot(np.arange(nTimesteps), objective[molecule] * np.ones(nTimesteps))

	ax.set_title("%s" % molecule, fontsize = 6)
	ax.set_ylabel("Counts", fontsize = 6)
	ax.set_xlabel("Timesteps", fontsize = 6)
	ax.spines['top'].set_visible(False)
	ax.spines['bottom'].set_visible(False)
	ax.xaxis.set_ticks_position('none')
	ax.tick_params(which = 'both', direction = 'out', labelsize = 6)
	ax.set_xticks([])

plt.subplots_adjust(hspace = 1, wspace = 1)
plotOutDir = current_path
plotOutFileName = "homeostatic_FBA_modular_FBA"

exportFigure(plt, plotOutDir, plotOutFileName, metadata = None)
plt.close()