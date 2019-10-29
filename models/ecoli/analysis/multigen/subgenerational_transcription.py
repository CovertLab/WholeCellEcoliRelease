"""
Generates plots related to the subgenerational transcription of genes.

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/12/2017
"""

from __future__ import absolute_import

import os
import cPickle

import numpy as np
import matplotlib.pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.utils.sparkline import whitePadSparklineAxis
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot

BUILD_CACHE = True
PLOT_GENES_OF_INTEREST = False
PLOT_DENOMINATOR_N_EACH_FREQ_GROUP = False
FIRST_N_GENS = 5

def remove_xaxis(axis):
	axis.spines["bottom"].set_visible(False)
	axis.tick_params(bottom = "off")
	axis.tick_params(axis = "x", labelbottom = 'off')
	axis.set_xlabel("")

class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(seedOutDir):
			raise Exception, "seedOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Get all cells
		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
		allDir = ap.get_cells()

		if len(allDir) <= 1:
			print "Skipping - this plot only runs for multigen sims"
			return

		sim_data = cPickle.load(open(simDataFile, "rb"))
		validation_data = cPickle.load(open(validationDataFile, "rb"))

		# Get mRNA data
		rnaIds = sim_data.process.transcription.rnaData["id"]
		isMRna = sim_data.process.transcription.rnaData["isMRna"]
		mRnaIndexes = np.where(isMRna)[0]
		mRnaIds = np.array([rnaIds[x] for x in mRnaIndexes])

		time = []
		time_eachGen = []
		transcribedBool = []
		simulatedSynthProbs = []
		transcriptionEvents = []

		for gen, simDir in enumerate(allDir):
			simOutDir = os.path.join(simDir, "simOut")

			if gen < FIRST_N_GENS:
				time += TableReader(os.path.join(simOutDir, "Main")).readColumn("time").tolist()
				time_eachGen.append(TableReader(os.path.join(simOutDir, "Main")).readColumn("time").tolist()[0])

			rnaSynthProb = TableReader(os.path.join(simOutDir, "RnaSynthProb"))
			simulatedSynthProb = np.mean(rnaSynthProb.readColumn("rnaSynthProb")[:, mRnaIndexes], axis = 0)
			rnaSynthProb.close()
			simulatedSynthProbs.append(simulatedSynthProb)

			bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
			moleculeIds = bulkMolecules.readAttribute("objectNames")
			mRnaIndexes_bulk = np.array([moleculeIds.index(x) for x in mRnaIds])
			moleculeCounts = bulkMolecules.readColumn("counts")[:, mRnaIndexes_bulk]
			bulkMolecules.close()
			moleculeCountsSumOverTime = moleculeCounts.sum(axis = 0)
			mRnasTranscribed = np.array([x != 0 for x in moleculeCountsSumOverTime])
			transcribedBool.append(mRnasTranscribed)

			rnapDataReader = TableReader(os.path.join(simOutDir, "RnapData"))
			rnaInitEvent = rnapDataReader.readColumn("rnaInitEvent")[:, mRnaIndexes]
			rnapDataReader.close()

			if gen == 0:
				transcriptionEvents = (rnaInitEvent != 0)
			elif gen < FIRST_N_GENS:
				transcriptionEvents = np.vstack((transcriptionEvents, (rnaInitEvent != 0)))
			else:
				pass

		time = np.array(time)
		time_eachGen.append(time[-1])
		time_eachGen = np.array(time_eachGen)
		transcribedBool = np.array(transcribedBool)
		simulatedSynthProbs = np.array(simulatedSynthProbs)

		indexingOrder = np.argsort(np.mean(simulatedSynthProbs, axis = 0))
		transcribedBoolOrdered = np.mean(transcribedBool, axis = 0)[indexingOrder]
		transcriptionEventsOrdered = transcriptionEvents[:, indexingOrder]
		mRnaIdsOrdered = mRnaIds[indexingOrder]

		alwaysPresentIndexes = np.where(transcribedBoolOrdered == 1.)[0]
		neverPresentIndexes = np.where(transcribedBoolOrdered == 0.)[0]
		sometimesPresentIndexes = np.array([x for x in np.arange(len(transcribedBoolOrdered)) if x not in alwaysPresentIndexes and x not in neverPresentIndexes])

		plt.style.use('seaborn-deep')
		color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
		color_always = color_cycle[2]
		color_never = color_cycle[4]
		color_subgen = color_cycle[0]
		colors = np.repeat(color_subgen, len(transcribedBoolOrdered))
		colors[alwaysPresentIndexes] = color_always
		colors[neverPresentIndexes] = color_never
		always = transcribedBoolOrdered[alwaysPresentIndexes]
		never = transcribedBoolOrdered[neverPresentIndexes]
		sometimes = transcribedBoolOrdered[sometimesPresentIndexes]

		alwaysTranscriptionEvents = []
		for i in alwaysPresentIndexes:
			v = (time[transcriptionEventsOrdered[:, i]] / 3600.).tolist()
			if transcriptionEventsOrdered[:, i].sum() == 0:
				v = [-1]
			alwaysTranscriptionEvents.append(v)

		neverTranscriptionEvents = []
		for i in neverPresentIndexes:
			v = (time[transcriptionEventsOrdered[:, i]] / 3600.).tolist()
			if transcriptionEventsOrdered[:, i].sum() == 0:
				v = [-1]
			neverTranscriptionEvents.append(v)

		sometimesTranscriptionEvents = []
		for i in sometimesPresentIndexes:
			v = (time[transcriptionEventsOrdered[:, i]] / 3600.).tolist()
			if transcriptionEventsOrdered[:, i].sum() == 0:
				v = [-1]
			sometimesTranscriptionEvents.append(v)

		# Plot
		plt.figure(figsize = (11, 8))
		scatterAxis = plt.subplot2grid((2, 4), (0, 0), colspan = 3, rowspan = 2)
		histAxis = plt.subplot2grid((2, 4), (0, 3), colspan = 1, rowspan = 2, sharey = scatterAxis)

		scatterAxis.scatter(np.arange(len(transcribedBoolOrdered)), transcribedBoolOrdered, marker = 'o', facecolors = colors, edgecolors = "none", s = 20)
		scatterAxis.set_xlim([0, len(transcribedBoolOrdered)])
		scatterAxis.set_ylim([0, 1])
		whitePadSparklineAxis(scatterAxis)

		N, bins, patches = histAxis.hist(transcribedBoolOrdered, bins = len(allDir) + 1, orientation = 'horizontal')
		for i in xrange(1, len(patches) - 1):
			plt.setp(patches[i], facecolor = "none", edgecolor = color_subgen)
		plt.setp(patches[0], facecolor = "none", edgecolor = color_never)
		plt.setp(patches[-1], facecolor = "none", edgecolor = color_always)
		whitePadSparklineAxis(histAxis)
		histAxis.xaxis.tick_bottom()
		histXmin, histXmax = histAxis.get_xlim()

		scatterAxis.set_ylim([-.01, 1.01])
		scatterAxis.set_yticklabels([])
		scatterAxis.set_xticklabels([])
		histAxis.set_xscale("log")
		histAxis.set_xlim([10, 2000])
		histAxis.set_yticklabels([])
		histAxis.set_xticklabels([])

		plt.subplots_adjust(wspace = 0.6, hspace = 0.4, right = 0.9, bottom = 0.1, left = 0.1, top = 0.9)
		exportFigure(plt, plotOutDir, plotOutFileName + "_frequency_histogram__clean", metadata)

		if PLOT_GENES_OF_INTEREST:
			## Identifying particular genes
			dcurId = "G7826_RNA[c]"
			clppId = "EG10158_RNA[c]"
			dcucId = "G6347_RNA[c]"
			dcurIndex = np.where(mRnaIdsOrdered == dcurId)[0]
			clppIndex = np.where(mRnaIdsOrdered == clppId)[0]
			dcucIndex = np.where(mRnaIdsOrdered == dcucId)[0]

			scatterAxis.scatter(dcurIndex, transcribedBoolOrdered[dcurIndex], facecolor = "orange", edgecolors = "none")
			scatterAxis.scatter(clppIndex, transcribedBoolOrdered[clppIndex], facecolor = "orange", edgecolors = "none")
			scatterAxis.scatter(dcucIndex, transcribedBoolOrdered[dcucIndex], facecolor = "orange", edgecolors = "none")
			scatterAxis.text(dcurIndex + 500, transcribedBoolOrdered[dcurIndex], "dcuR", fontsize = 18)
			scatterAxis.text(clppIndex + 500, transcribedBoolOrdered[clppIndex], "clpP", fontsize = 18)
			scatterAxis.text(dcucIndex + 500, transcribedBoolOrdered[dcucIndex], "dcuC", fontsize = 18)
			exportFigure(plt, plotOutDir, plotOutFileName + "_frequency_histogram__clean__genes", metadata)

		plt.suptitle("Frequency of observing at least 1 transcript per generation", fontsize = 14)
		scatterAxis.set_xlabel("Genes ordered by simulated synthesis probability", fontsize = 12)
		histAxis.text(histXmax * 1.6, 0, "%s genes\n(%0.1f%%)" % (len(never), 100. * (len(never) / float(len(transcribedBoolOrdered)))), fontsize = 14, verticalalignment = "center", color = color_never)
		histAxis.text(histXmax * 1.6, 1, "%s genes\n(%0.1f%%)" % (len(always), 100. * (len(always) / float(len(transcribedBoolOrdered)))), fontsize = 14, verticalalignment = "center", color = color_always)
		histAxis.text(histXmax * 1.6, 0.5, "%s genes\n(%0.1f%%)" % (len(sometimes), 100. * (len(sometimes) / float(len(transcribedBoolOrdered)))), fontsize = 14, verticalalignment = "center", color = color_subgen)
		scatterAxis.set_yticklabels([0, 1])
		scatterAxis.set_xticklabels([0, len(transcribedBoolOrdered)])
		histAxis.set_xticks([10, 2000])
		histAxis.set_xticklabels([10, 2000])
		histAxis.set_yticklabels([0, 1])
		exportFigure(plt, plotOutDir, plotOutFileName + "_frequency_histogram", metadata)
		plt.close("all")

		plt.figure(figsize = (16, 8))
		alwaysAxis = plt.subplot(2, 1, 1)
		sometimesAxis = plt.subplot(2, 1, 2, sharex = alwaysAxis)

		alwaysAxis.eventplot(alwaysTranscriptionEvents, orientation = "horizontal", linewidths = 2., linelengths = 4., color = color_always)
		alwaysAxis.set_xlim([0, time[-1] / 3600.])
		alwaysAxis.set_ylim([-1, len(always)])
		alwaysAxis.set_xticks([])
		alwaysAxis.set_yticks([])
		alwaysAxis.tick_params(top = "off")
		alwaysAxis.tick_params(bottom = "off")
		alwaysAxis.tick_params(axis = "x", labelbottom = 'off')
		
		sometimesAxis.eventplot(sometimesTranscriptionEvents, orientation = "horizontal", linewidths = 2., linelengths = 4., color = color_subgen)
		sometimesAxis.set_xlim([0, time[-1] / 3600.])
		sometimesAxis.set_ylim([-1, len(sometimes)])
		sometimesAxis.set_xticks([0, time[-1] / 3600.])
		sometimesAxis.set_yticks([])
		sometimesAxis.tick_params(top = "off")
		sometimesAxis.tick_params(which = 'both', direction = 'out', labelsize = 12)
		time_eachGen = np.array(time_eachGen)
		sometimesAxis.set_xticks(time_eachGen / 3600.)
		sometimesAxis.set_xticklabels([])

		plt.subplots_adjust(wspace = 0, hspace = 0, right = 0.9, bottom = 0.1, left = 0.1, top = 0.9)
		exportFigure(plt, plotOutDir, plotOutFileName + "_eventplot__clean", metadata)

		plt.suptitle("Transcription initiation events", fontsize = 14)
		alwaysAxis.set_ylabel("Freq. = 1", fontsize = 14)
		sometimesAxis.set_ylabel("0 < Freq. < 1", fontsize = 14)
		sometimesAxis.set_xlabel("Time (gens)", fontsize = 14)
		sometimesAxis.set_xticklabels(np.arange(FIRST_N_GENS + 1))
		exportFigure(plt, plotOutDir, plotOutFileName + "_eventplot", metadata)
		plt.close("all")

		# Essential genes figure
		nRed = len(never)
		nGreen = len(sometimes)
		nBlue = len(always)

		plotRed = 0
		plotGreen = 0
		plotBlue = 0

		essentialGenes_rna = validation_data.essentialGenes.essentialRnas
		for g in essentialGenes_rna:
			i = np.where(mRnaIdsOrdered == str(g))[0][0]
			f = transcribedBoolOrdered[i]

			if f == 0.0:
				plotRed += 1
			elif f == 1.0:
				plotBlue += 1
			else:
				plotGreen += 1
		xloc = np.arange(3)
		width = 0.8

		plt.figure()
		ax = plt.subplot(1, 1, 1)
		ax.bar(xloc + width, [plotRed / float(len(essentialGenes_rna)), plotGreen / float(len(essentialGenes_rna)), plotBlue / float(len(essentialGenes_rna))], width, color = [color_never, color_subgen, color_always], edgecolor ="none")
		whitePadSparklineAxis(ax)
		ax.spines["left"].set_position(("outward", 0))
		ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
		ax.set_yticklabels(["0%", "20%", "40%", "60%", "80%"])
		ax.set_ylabel("Percentage of essential genes")
		ax.set_xticks(xloc + 1.5 * width)
		ax.set_xticklabels(["%s / %s" % (plotRed, len(essentialGenes_rna)), "%s / %s" % (plotGreen, len(essentialGenes_rna)), "%s / %s" % (plotBlue, len(essentialGenes_rna))])
		ax.set_xlabel("Total number of essential genes: %s" % len(essentialGenes_rna))
		plt.subplots_adjust(right = 0.9, bottom = 0.15, left = 0.2, top = 0.9)
		exportFigure(plt, plotOutDir, plotOutFileName + "_essential_genes", metadata)

		ax.set_yticklabels([])
		ax.set_ylabel("")
		remove_xaxis(ax)
		exportFigure(plt, plotOutDir, plotOutFileName + "_essential_genes__clean", metadata)

		if PLOT_DENOMINATOR_N_EACH_FREQ_GROUP:	
			plt.figure()
			ax = plt.subplot(1, 1, 1)
			ax.bar(xloc + width, [plotRed / float(nRed), plotGreen / float(nGreen), plotBlue / float(nBlue)], width, color = [color_never, color_subgen, color_always], edgecolor ="none")
			whitePadSparklineAxis(ax)
			ax.spines["left"].set_position(("outward", 0))
			ax.set_yticks([0.0, 0.1, 0.2])
			ax.set_yticklabels(["0%", "10%", "20%"])
			ax.set_ylabel("Percentage of genes that are essential genes")
			ax.set_xticks(xloc + 1.5 * width)
			ax.set_xticklabels(["%s / %s" % (plotRed, nRed), "%s / %s" % (plotGreen, nGreen), "%s / %s" % (plotBlue, nBlue)])
			plt.subplots_adjust(right = 0.9, bottom = 0.1, left = 0.2, top = 0.9)
			exportFigure(plt, plotOutDir, plotOutFileName + "_essential_genes_v2", metadata)
		plt.close()

		# Gene annotation/antibiotics figure
		geneFunctions = validation_data.geneFunctions.geneFunctions
		unknown = {"r": 0, "g": 0, "b": 0}
		resistance = {"r": 0, "g": 0, "b": 0}
		for frameID, function_ in geneFunctions.iteritems():
			if function_ in ["Unknown function", "Unclear/under-characterized"]:
				i = np.where([frameID in x for x in mRnaIdsOrdered])[0][0]
				f = transcribedBoolOrdered[i]
				if f == 0.0:
					unknown["r"] += 1
				elif f == 1.0:
					unknown["b"] += 1
				else:
					unknown["g"] += 1

			elif function_ in ["Antibiotic resistance", "Toxin/antitoxin"]:
				i = np.where([frameID in x for x in mRnaIdsOrdered])[0][0]
				f = transcribedBoolOrdered[i]
				if f == 0.0:
					resistance["r"] += 1
				elif f == 1.0:
					resistance["b"] += 1
				else:
					resistance["g"] += 1

		nUnknown = np.sum([unknown[x] for x in ["r", "g", "b"]])
		plt.figure()
		ax = plt.subplot(1, 1, 1)
		ax.bar(xloc + width, [unknown["r"] / float(nUnknown), unknown["g"] / float(nUnknown), unknown["b"] / float(nUnknown)], width, color = [color_never, color_subgen, color_always], edgecolor ="none")
		whitePadSparklineAxis(ax)
		ax.spines["left"].set_position(("outward", 0))
		ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
		ax.set_yticklabels(["0%", "20%", "40%", "60%", "80%"])
		ax.set_ylabel("Percentage of poorly understood / annotated genes")
		ax.set_xticks(xloc + 1.5 * width)
		ax.set_xticklabels(["%s / %s" % (unknown["r"], nUnknown), "%s / %s" % (unknown["g"], nUnknown), "%s / %s" % (unknown["b"], nUnknown)])
		ax.set_xlabel("Total number of poorly understood / annotated genes: %s" % nUnknown)
		plt.subplots_adjust(right = 0.9, bottom = 0.15, left = 0.2, top = 0.9)
		exportFigure(plt, plotOutDir, plotOutFileName + "_unannotated", metadata)

		ax.set_yticklabels([])
		ax.set_ylabel("")
		remove_xaxis(ax)
		exportFigure(plt, plotOutDir, plotOutFileName + "_unannotated__clean", metadata)

		if PLOT_DENOMINATOR_N_EACH_FREQ_GROUP:
			plt.figure()
			ax = plt.subplot(1, 1, 1)
			ax.bar(xloc + width, [unknown["r"] / float(nRed), unknown["g"] / float(nGreen), unknown["b"] / float(nBlue)], width, color = [color_never, color_subgen, color_always], edgecolor ="none")
			whitePadSparklineAxis(ax)
			ax.spines["left"].set_position(("outward", 0))
			ax.set_yticks([0.0, 0.2, 0.4, 0.6])
			ax.set_yticklabels(["0%", "20%", "40%", "60%"])
			ax.set_ylabel("Percentage of genes that are poorly understood / annotated")
			ax.set_xticks(xloc + 1.5 * width)
			ax.set_xticklabels(["%s / %s" % (unknown["r"], nRed), "%s / %s" % (unknown["g"], nGreen), "%s / %s" % (unknown["b"], nBlue)])
			plt.subplots_adjust(right = 0.9, bottom = 0.1, left = 0.2, top = 0.9)
			exportFigure(plt, plotOutDir, plotOutFileName + "_v2", metadata)
		plt.close()


		nResistance = float(np.sum([resistance[x] for x in ["r", "g", "b"]]))
		plt.figure()
		ax = plt.subplot(1, 1, 1)
		ax.bar(xloc + width, [resistance["r"] / nResistance, resistance["g"] / nResistance, resistance["b"] / nResistance], width, color = [color_never, color_subgen, color_always], edgecolor ="none")
		whitePadSparklineAxis(ax)
		ax.spines["left"].set_position(("outward", 0))
		ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
		ax.set_yticklabels(["0%", "20%", "40%", "60%", "80%"])
		ax.set_ylabel("Percentage of antibiotic related genes")
		ax.set_xticks(xloc + 1.5 * width)
		ax.set_xticklabels(["%s / %s" % (resistance["r"], int(nResistance)), "%s / %s" % (resistance["g"], int(nResistance)), "%s / %s" % (resistance["b"], int(nResistance))])
		ax.set_xlabel("Total number of antibiotic related genes: %s" % int(nResistance))
		plt.subplots_adjust(right = 0.9, bottom = 0.15, left = 0.2, top = 0.9)
		exportFigure(plt, plotOutDir, plotOutFileName + "_antibiotic", metadata)

		ax.set_yticklabels([])
		ax.set_ylabel("")
		remove_xaxis(ax)
		exportFigure(plt, plotOutDir, plotOutFileName + "_antibiotic__clean", metadata)

		if PLOT_DENOMINATOR_N_EACH_FREQ_GROUP:
			plt.figure()
			ax = plt.subplot(1, 1, 1)
			ax.bar(xloc + width, [resistance["r"] / float(nRed), resistance["g"] / float(nGreen), resistance["b"] / float(nBlue)], width, color = [color_never, color_subgen, color_always], edgecolor ="none")
			whitePadSparklineAxis(ax)
			ax.spines["left"].set_position(("outward", 0))
			ax.set_yticks([0.0, 0.025])
			ax.set_yticklabels(["0%", "2.5%"])
			ax.set_ylabel("Percentage of genes that are antibiotic related")
			ax.set_xticks(xloc + 1.5 * width)
			ax.set_xticklabels(["%s / %s" % (resistance["r"], nRed), "%s / %s" % (resistance["g"], nGreen), "%s / %s" % (resistance["b"], nBlue)])
			plt.subplots_adjust(right = 0.9, bottom = 0.1, left = 0.2, top = 0.9)
			exportFigure(plt, plotOutDir, plotOutFileName + "_antibiotic_v2", metadata)
		plt.close()


if __name__ == "__main__":
	Plot().cli()
