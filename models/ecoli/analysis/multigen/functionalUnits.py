"""
Groups subgenerationally transcribed genes by:
	- monomeric vs complexed functional units
	- always present vs 0 at least once
	- essential vs non-essential

Notes:
	- When determining the presence of functional units, the first 5 timesteps
	of the 1st generation are skipped because complexes are not yet formed
	during these steps. See SKIP_TIMESTEPS.
	- This script can output a summary of monomer and complex IDs found to be
	subgenerationally expressed. Table S4 was generated from this output.
	See WRITE_TO_FILE.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 7/12/2018
"""

from __future__ import absolute_import

import os
import cPickle
import numpy as np
import matplotlib.pyplot as plt
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot
from wholecell.utils.sparkline import whitePadSparklineAxis

WRITE_TO_FILE = True
SKIP_TIMESTEPS = 5

def plot(proteinIds, zeroAtLeastOnce, essential, ax, axEssential, xloc, width, functionalUnitIds, highlight_color):
	if not len(functionalUnitIds):
		notAlwaysPresent = zeroAtLeastOnce.sum()
		alwaysPresent = len(proteinIds) - notAlwaysPresent

		# Count number of essential genes
		notAlwaysPresentEssential = np.logical_and(zeroAtLeastOnce, essential).sum()
		alwaysPresentEssential = np.logical_and(np.logical_not(zeroAtLeastOnce), essential).sum()

	else:
		functionalUnitIndices = [proteinIds.index(x) for x in functionalUnitIds]
		notAlwaysPresent = zeroAtLeastOnce[functionalUnitIndices].sum()
		alwaysPresent = len(functionalUnitIds) - notAlwaysPresent

		# Count number of essential genes
		notAlwaysPresentEssential = np.logical_and(zeroAtLeastOnce[functionalUnitIndices], essential[functionalUnitIndices]).sum()
		alwaysPresentEssential = np.logical_and(np.logical_not(zeroAtLeastOnce[functionalUnitIndices]), essential[functionalUnitIndices]).sum()

	# Plot
	bar = ax.bar(xloc + width, [alwaysPresent, notAlwaysPresent], width, color = "lightgray")
	axEssential.bar(xloc + width, [alwaysPresent, notAlwaysPresent], width, color = "lightgray")
	barEssential = axEssential.bar(xloc + width, [alwaysPresentEssential, notAlwaysPresentEssential], width, color = highlight_color)

	# Format
	for axis in [ax, axEssential]:
		whitePadSparklineAxis(axis)
		axis.set_xticklabels(["Always\npresent", "0 at\nleast once"])
		axis.set_xticks(xloc + width)
	ax.set_yticks([x.get_height() for x in bar])
	axRight = axEssential.twinx()
	axRight.set_yticks([x[1].get_height() for x in [bar, barEssential]])
	axEssential.set_yticks([x[0].get_height() for x in [bar, barEssential]])
	axRight.set_ylim(ax.get_ylim())
	axRight.spines["left"].set_visible(False)
	axRight.spines["top"].set_visible(False)
	return

class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(seedOutDir):
			raise Exception, "seedOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		## Identify sub-genenerationally transcribed genes
		# Get all cells
		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
		allDir = ap.get_cells()

		# Load sim data
		sim_data = cPickle.load(open(simDataFile, "rb"))
		rnaIds = sim_data.process.transcription.rnaData["id"]
		geneIds = sim_data.process.transcription.rnaData["geneId"]

		# For each generation
		nonzeroSumRnaCounts_allGens = []
		for i, simDir in enumerate(allDir):
			simOutDir = os.path.join(simDir, "simOut")

			# Read counts of transcripts
			bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
			if i == 0:
				moleculeIds = bulkMolecules.readAttribute("objectNames")
				rnaIndices = np.array([moleculeIds.index(x) for x in rnaIds])
			rnaCounts = bulkMolecules.readColumn("counts")[:, rnaIndices]
			bulkMolecules.close()

			# Sum counts over timesteps
			sumRnaCounts = rnaCounts.sum(axis = 0)

			# Flag where the sum is nonzero (True if nonzero, False if zero)
			nonzeroSumRnaCounts = sumRnaCounts != 0
			nonzeroSumRnaCounts_allGens.append(nonzeroSumRnaCounts)

		# Average (mean) over generations
		nonzeroSumRnaCounts_allGens = np.array(nonzeroSumRnaCounts_allGens)
		avgRnaCounts = nonzeroSumRnaCounts_allGens.mean(axis = 0)

		# Identify subgenerationally transcribed genes
		subgenRnaIndices = np.where(np.logical_and(avgRnaCounts != 0., avgRnaCounts != 1.))[0]
		subgenRnaIds = rnaIds[subgenRnaIndices]
		subgenMonomerIndices = [np.where(sim_data.process.translation.monomerData["rnaId"] == x)[0][0] for x in subgenRnaIds]
		subgenMonomerIds = sim_data.process.translation.monomerData["id"][subgenMonomerIndices]

		# Identify subgenerationally transcribed genes that function as monomers
		complexationMoleculeNames = sim_data.process.complexation.moleculeNames
		subgenMonomerOnlyIds = [x for x in subgenMonomerIds if x not in complexationMoleculeNames]
		subgenMonomersInComplexesIds = [x for x in subgenMonomerIds if x in complexationMoleculeNames]

		## Identify complexes that subgenerationally transcribed genes participate in
		subgenComplexIds = []
		for complexId in sim_data.process.complexation.complexNames:
			subunitIds = sim_data.process.complexation.getMonomers(complexId)["subunitIds"]
			if np.any([x in subgenMonomerIds for x in subunitIds]):
				subgenComplexIds.append(complexId)
		subgenComplexIds = list(set(subgenComplexIds))

		## Identify functional units that have a 0 count for at least one timestep
		proteinIds = list(set(subgenMonomerOnlyIds + subgenComplexIds))
		if not len(proteinIds):
			print "Returned -- No subgenerational functional units were found."
			return
		zeroAtLeastOnce = np.zeros(len(proteinIds), dtype = bool)
		for i, simDir in enumerate(allDir):
			simOutDir = os.path.join(simDir, "simOut")

			# Read counts of proteins
			bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
			if i == 0:
				moleculeIds = bulkMolecules.readAttribute("objectNames")
				proteinIndices = np.array([moleculeIds.index(x) for x in proteinIds])
				proteinCounts = bulkMolecules.readColumn("counts")[SKIP_TIMESTEPS:, proteinIndices]
			else:
				proteinCounts = bulkMolecules.readColumn("counts")[:, proteinIndices]
			bulkMolecules.close()

			# Flag proteins with a minimum count of 0
			minProteinCounts = np.min(proteinCounts, axis = 0)
			zeroAtLeastOnce[minProteinCounts == 0] = True

		## Identify essential functional units
		# Load validation data
		validation_data = cPickle.load(open(validationDataFile, "rb"))
		essentialGenes_genes = validation_data.essentialGenes.essentialGenes
		essentialGenes_monomers = validation_data.essentialGenes.essentialProteins

		# Flag essential monomers
		essential = np.zeros(len(proteinIds), dtype = bool)
		essential[[x in essentialGenes_monomers for x in proteinIds]] = True

		# Flag essential complexes
		subgenMonomerIdToComplexIds = {}
		for complexId in subgenComplexIds:
			subunitIds = sim_data.process.complexation.getMonomers(complexId)["subunitIds"]
			isEssential = np.any([(x in essentialGenes_monomers) and (x in subgenMonomerIds) for x in subunitIds])

			if isEssential:
				complexIndex = proteinIds.index(complexId)
				essential[complexIndex] = True

			# Add unadded complexes to subgenMonomerIdToComplexIds
			for monomerId in subunitIds:
				if monomerId not in subgenMonomerIdToComplexIds.keys():
					subgenMonomerIdToComplexIds[monomerId] = []
				subgenMonomerIdToComplexIds[monomerId].append(complexId)

		## Plot
		nrows = 2
		ncols = 3
		fig, axesList = plt.subplots(nrows, ncols, figsize = (11, 8.5))
		[[axMonomers, axComplexes, axTotal], [axMonomersEssential, axComplexesEssential, axTotalEssential]] = axesList
		xloc = np.arange(2)
		width = 0.75
		plt.style.use('seaborn-deep')
		highlight_color = plt.rcParams['axes.prop_cycle'].by_key()['color'][0]

		# Plot subgenerational genes that don't form complexes
		args = {"proteinIds": proteinIds,
			"zeroAtLeastOnce": zeroAtLeastOnce,
			"essential": essential,
			"ax": axMonomers,
			"axEssential": axMonomersEssential,
			"xloc": xloc,
			"width": width,
			"functionalUnitIds": subgenMonomerOnlyIds,
			"highlight_color": highlight_color,
			}
		plot(**args)
		axMonomers.set_title("{0} monomeric\nfunctional units".format(len(subgenMonomerOnlyIds)))
		axMonomers.set_ylabel("# functional units")
		axMonomersEssential.set_ylabel("# functional units\n(essential genes in blue)")

		# Plot subgenerational genes that form complexes
		args["ax"] = axComplexes
		args["axEssential"] = axComplexesEssential
		args["functionalUnitIds"] = subgenComplexIds
		plot(**args)
		axComplexes.set_title("{0} complexed\nfunctional units".format(len(subgenComplexIds)))

		# Plot subgenenerational functional units
		args["ax"] = axTotal
		args["axEssential"] = axTotalEssential
		args["functionalUnitIds"] = []
		plot(**args)
		axTotal.set_title("{0} + {1} = {2}\nfunctional units".format(len(subgenMonomerOnlyIds), len(subgenComplexIds), len(proteinIds)))

		plt.subplots_adjust(wspace = 0.3)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")

		# Plot unlabeled version of final panel for Fig. 4E
		fig, ax = plt.subplots(1, 1, figsize = (5, 10))
		notAlwaysPresent = zeroAtLeastOnce.sum()
		alwaysPresent = len(proteinIds) - notAlwaysPresent
		notAlwaysPresentEssential = np.logical_and(zeroAtLeastOnce, essential).sum()
		alwaysPresentEssential = np.logical_and(np.logical_not(zeroAtLeastOnce), essential).sum()
		barAll = ax.bar(xloc + width, [alwaysPresent, notAlwaysPresent], width, color = "lightgray")
		barEssential = ax.bar(xloc + width, [alwaysPresentEssential, notAlwaysPresentEssential], width, color = highlight_color)

		# Format
		whitePadSparklineAxis(ax)
		ax.set_xticks(xloc + width)
		ax.set_xticklabels([])
		ax.set_yticks([0] + [x[0].get_height() for x in [barAll, barEssential]])

		axRight = ax.twinx()
		axRight.set_ylim(ax.get_ylim())
		axRight.set_yticks([0] + [x[1].get_height() for x in [barAll, barEssential]])

		axRight.spines["left"].set_visible(False)
		axRight.spines["top"].set_visible(False)
		axRight.spines["right"].set_position(("outward", 10))
		axRight.spines["bottom"].set_position(("outward", 10))
		for axis in [ax, axRight]:
			axis.set_yticklabels([])
			axis.tick_params(length = 10.0)

		plt.subplots_adjust(bottom = 0.3, top = 0.7, left = 0.2, right = 0.8)
		exportFigure(plt, plotOutDir, "{0}_clean".format(plotOutFileName), metadata)
		plt.close("all")

		if WRITE_TO_FILE:
			with open(os.path.join(plotOutDir, "%s.tsv" % plotOutFileName), "wb") as f:
				f.write("gene\trna\tmonomer\tcomplex(es)\tzero_at_least_once\tessential\n")
				for rnaId, monomerId in zip(subgenRnaIds, subgenMonomerIds):
					# Get gene Id
					rnaIndex = np.where(rnaIds == rnaId)[0][0]
					geneId = geneIds[rnaIndex]

					# Get IDs of complex(es)
					if monomerId in subgenMonomerIdToComplexIds:
						complexIds = subgenMonomerIdToComplexIds[monomerId]
						functionalUnitIds = complexIds
					else:
						complexIds = None
						functionalUnitIds = [monomerId]

					# Were funcational units 0 at least once
					isZero = []
					for functionalUnitId in functionalUnitIds:
						functionalUnitIndex = proteinIds.index(functionalUnitId)
						isZero.append(zeroAtLeastOnce[functionalUnitIndex])

					# Is this gene essential
					isEssential = geneId in essentialGenes_genes

					# Write
					f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(geneId, rnaId, monomerId, complexIds, isZero, isEssential))


if __name__ == "__main__":
	Plot().cli()
