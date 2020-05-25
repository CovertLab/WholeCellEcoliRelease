"""
@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 11/13/2017
"""

from __future__ import absolute_import, division, print_function

import os
import cPickle

import numpy as np
import matplotlib.pyplot as plt

from wholecell.utils.sparkline import whitePadSparklineAxis
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Check if cache from rnaVsProteinPerCell.py exists
		if os.path.exists(os.path.join(plotOutDir, "rnaVsProteinPerCell_alltimesteps.cPickle")):
			rnaVsProteinPerCell = cPickle.load(open(os.path.join(plotOutDir, "rnaVsProteinPerCell_alltimesteps.cPickle"), "rb"))
			avgProteinCounts_forAllCells = rnaVsProteinPerCell["protein"]
			avgProteinCounts_perCell = avgProteinCounts_forAllCells / float(32)
		else:
			print("Requires rnaVsProteinPerCell.cPickle from rnaVsProteinPerCell.py")
			return

		# Check if cache from figure5B_E_F_G.py exist
		if os.path.exists(os.path.join(plotOutDir, "figure5B.pickle")):
			figure5B_data = cPickle.load(open(os.path.join(plotOutDir, "figure5B.pickle"), "rb"))
			colors = figure5B_data["colors"]
			mrnaIds = figure5B_data["id"].tolist()
		else:
			print("Requires figure5B.pickle from figure5B_E_F_G.py")
			return

		# Check if cache functionalUnits.py exist
		if os.path.exists(os.path.join(plotOutDir, "functionalUnits.cPickle")):
			functionalUnits_data = cPickle.load(open(os.path.join(plotOutDir, "functionalUnits.cPickle"), "rb"))
			minProteinCounts = functionalUnits_data["minProteinCounts"]
			monomersInManyComplexes_dict = functionalUnits_data["monomersInvolvedInManyComplexes_dict"]
		else:
			print("Requires functionalUnits.cPickle from functionalUnits.py")
			return

		# Check if cache ratioFinalToInitialCountMultigen.pickle from figure5_c.py (cohort analysis) exists
		cohortPlotOutDir = "out/SET_A/20170407.100507.741102__SET_A_32_gens_8_seeds_basal_with_growth_noise_and_D_period/wildtype_000000/plotOut/" # todo, generalize
		if os.path.exists(os.path.join(cohortPlotOutDir, "ratioFinalToInitialCountMultigen.pickle")):
			ratioFinalToInitialCountMultigen = cPickle.load(open(os.path.join(plotOutDir,"ratioFinalToInitialCountMultigen.pickle"), "rb"))
		else:
			print("Requires ratioFinalToInitialCountMultigen.pickle from figure5_c.py")
			return

		# Load sim data
		sim_data = cPickle.load(open(simDataFile, "rb"))
		rnaIds = sim_data.process.transcription.rnaData["id"][sim_data.relation.rnaIndexToMonomerMapping] # orders rna IDs to match monomer IDs
		rnaIds = rnaIds.tolist()
		ids_complexation = sim_data.process.complexation.moleculeNames # Complexe of proteins, and protein monomers
		ids_complexation_complexes = sim_data.process.complexation.ids_complexes # Only complexes
		ids_translation = sim_data.process.translation.monomerData["id"].tolist() # Only protein monomers

		# ID subgenerational rnas
		subgenerational_indices = np.where(colors == "b")[0]
		subgenerational_rnaIds = np.array(mrnaIds)[subgenerational_indices]

		# Get min count of monomers for each subgenerational rna
		rnaIndices = [rnaIds.index(x) for x in subgenerational_rnaIds]
		minProteinCounts = minProteinCounts[rnaIndices]

		# Get min count of functional units (ie. monomers involved in complexes) for each subgeneration rna
		subgenerational_monomerIds = np.array(ids_translation)[rnaIndices]
		subgenerational_monomersInComplex = [x for x in subgenerational_monomerIds if x in monomersInManyComplexes_dict]
		subgenerational_monomersNotInComplex = [x for x in subgenerational_monomerIds if x not in subgenerational_monomersInComplex]
		complexMinCounts_dict = {}
		for i in monomersInManyComplexes_dict.keys():
			if i not in subgenerational_monomerIds:
				continue
			for j in monomersInManyComplexes_dict[i].keys():
				count = monomersInManyComplexes_dict[i][j]
				complexMinCounts_dict[j] = count

		# Check whether all subunits of subgenerational complexes are subgenerational monomers
		nComplexesWNonSubGSubunits = 0
		for i in complexMinCounts_dict.keys():
			subunitIds = sim_data.process.complexation.getMonomers(i)["subunitIds"]
			if not np.all([x in subgenerational_monomerIds for x in subunitIds]):
				nComplexesWNonSubGSubunits += 1

		# Get min count of each monomer
		noncomplex_zero = []
		noncomplex_one = []
		noncomplex_zero_monomerIDs = []
		noncomplex_one_monomerIDs = []
		for i, minProteinCount in enumerate(minProteinCounts):
			if minProteinCount == np.inf:
				continue
			if minProteinCount == 0:
				noncomplex_zero.append(i)
				noncomplex_zero_monomerIDs.append(subgenerational_monomerIds[i])
			else:
				noncomplex_one.append(i)
				noncomplex_one_monomerIDs.append(subgenerational_monomerIds[i])
		complex_zero = []
		complex_one = []
		for i in complexMinCounts_dict.keys():
			count = complexMinCounts_dict[i]
			if count == 0:
				complex_zero.append(i)
			else:
				complex_one.append(i)

		# Get indices of subgenerational genes (in order of ids_translation)
		subG_indices = [ids_translation.index(x) for x in subgenerational_monomerIds]
		subG_noncomplex_indices = [ids_translation.index(x) for x in subgenerational_monomersNotInComplex]
		subG_noncomplex_zero_indices = [ids_translation.index(x) for x in noncomplex_zero_monomerIDs]
		subG_noncomplex_one_indices = [ids_translation.index(x) for x in noncomplex_one_monomerIDs]

		# Get average fold change
		avgFC = np.average(ratioFinalToInitialCountMultigen[:, subgenerational_indices], axis = 0)
		avgFC_zero = np.average(ratioFinalToInitialCountMultigen[:, subG_noncomplex_zero_indices], axis = 0)
		avgFC_one = np.average(ratioFinalToInitialCountMultigen[:, subG_noncomplex_one_indices], axis = 0)

		# Plot
		plt.figure(figsize = (12, 12))
		nrows = 3
		ncols = 3
		ax = plt.subplot2grid((nrows, ncols), (0, 0))
		ax1 = plt.subplot2grid((nrows, ncols), (0, 1))
		ax2 = plt.subplot2grid((nrows, ncols), (0, 2))
		ax3 = plt.subplot2grid((nrows, ncols), (1, 0), colspan = ncols - 1)
		ax4 = plt.subplot2grid((nrows, ncols), (2, 0), colspan = ncols - 1, sharey = ax3)
		ax5 = plt.subplot2grid((nrows, ncols), (1, 2))
		ax6 = plt.subplot2grid((nrows, ncols), (2, 2))
		axesList = [ax, ax1, ax2, ax3, ax4, ax5, ax6]
		barAxesList = [ax, ax1, ax2]

		# Plot text
		nSubG = len(subgenerational_monomerIds)
		nSubGMonomers = len(subG_noncomplex_indices)
		nSubGMonomersInComplexes = nSubG - nSubGMonomers
		nSubGComplexes = len(complexMinCounts_dict.keys())
		nSubGFunctionalUnits = nSubGMonomers + nSubGComplexes

		text = "%s total subgenerational genes" % nSubG
		text += "\n-- (panel A) %s function as monomers" % nSubGMonomers
		text += "\n-- (panel B) %s are subunits to %s unique complexes" % (nSubGMonomersInComplexes, nSubGComplexes)
		text += "\n-- (panel C) sum of panels A and B"
		text += "\n-- (panel D) and (panel E):"
		text += "\n        -- colors are preserved from panel A"
		text += "\n        -- uses FC data from Figure 6D (takes an avg over the gens so that there is 1 data point per gene)"
		text += "\n        -- only represents monomers (ie. genes from panel A)"
		text += "\n           (because Figure 6D computes FCs after dissociating complexes into monomers)"
		text += "\n-- (panel F) is a linear y-scale version of panel D"
		text += "\n-- (panel G) is a linear y-scale version of panel E"
		plt.suptitle(text, multialignment = "left", fontsize = 10)

		# Plot bar graphs
		xloc = np.arange(2)
		width = 0.75

		bars = ax.bar(xloc + width, [len(noncomplex_one), len(noncomplex_zero)], width, edgecolor = "none")
		bars1 = ax1.bar(xloc + width, [len(complex_one), len(complex_zero)], width, edgecolor = "none")
		bars2 = ax2.bar(xloc + width, [len(noncomplex_one) + len(complex_one), len(noncomplex_zero) + len(complex_zero)], width, edgecolor = "none")
		ax.set_ylabel("# monomers")
		ax1.set_ylabel("# complexes")
		ax2.set_ylabel("# functional units")
		ax.set_title("%s monomers" % nSubGMonomers, fontsize = 10, y = 1.08)
		ax1.set_title("%s complexes" % nSubGComplexes, fontsize = 10, y = 1.08)
		ax2.set_title("%s + %s = %s functional units" % (nSubGMonomers, nSubGComplexes, nSubGFunctionalUnits), fontsize = 10, y = 1.08)
		for ax, bar in zip(barAxesList, [bars, bars1, bars2]):
			bar[0].set_alpha(0.5)
			bar[0].set_edgecolor("none")
			whitePadSparklineAxis(ax)
			ax.spines["left"].set_position(("outward", 0))
			ax.set_yticks([bar[0].get_height(), bar[1].get_height()])
			ax.set_xticks(xloc + 1.5 * width)
			ax.set_xticklabels(["always\npresent", "0 at\nleast once"])

		# Plot fold change histogram
		binRange = np.arange(0, 25, 0.4)
		n3, bins3, patches3 = ax3.hist(avgFC_zero, bins = binRange, log = True)
		n4, bins4, patches4 = ax4.hist(avgFC_one, alpha = 0.5, bins = binRange, log = True)
		ax5.hist(avgFC_zero, bins = binRange)
		ax6.hist(avgFC_one, bins = binRange, alpha = 0.5)
		xmin = min(np.floor(avgFC_zero.min()), np.floor(avgFC_one.min()))
		xmax = max(np.ceil(avgFC_zero.max()), np.ceil(avgFC_one.max()))
		for ax, bin_, n in zip([ax3, ax4, ax5, ax6], [bins3, bins4, bins3, bins4], [n3, n4, n3, n4]):
			ax.set_xlim([xmin, xmax])
			ax.set_xlabel("<Fold Change>")
			largestBin = np.argsort(n)[-1]
			ax.set_xticks([ax.get_xlim()[0], np.average(bin_[largestBin: largestBin + 2]), ax.get_xlim()[1]])

		for ax, label in zip(axesList, ['A', 'B', 'C', 'D', "E", "F", "G"]):
			ax.text(-0.1, 1.35, label, transform=ax.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
		plt.subplots_adjust(hspace = 0.6, wspace = 0.6, left = 0.1, bottom = 0.1, top = 0.75, right = 0.9)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")

		# Plot clean version of panel C for paper
		fig, ax = plt.subplots(1, 1, figsize = (5, 5))
		bars2 = ax.bar(xloc + width, [len(noncomplex_one) + len(complex_one), len(noncomplex_zero) + len(complex_zero)], width, edgecolor = "none")
		bars2[0].set_alpha(0.5)
		bars2[0].set_edgecolor("none")
		whitePadSparklineAxis(ax)
		ax.spines["left"].set_position(("outward", 0))
		ax.set_yticks([bars2[0].get_height(), bars2[1].get_height()])
		ax.set_xticks(xloc + 1.5 * width)
		ax.set_xticklabels([])
		ax.set_yticklabels([])
		exportFigure(plt, plotOutDir, plotOutFileName + "_C", metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
