from __future__ import absolute_import
from __future__ import division

import os

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as patches
import cPickle

from wholecell.io.tablereader import TableReader
from wholecell.utils import units
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot

class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile,
		validationDataFile, metadata):
		if not os.path.isdir(seedOutDir):
			raise Exception, "seedOutDir does not currently \
			exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)

		allDirs = ap.get_cells()

		###################
		# Initialize Figure:
		num_subplots = 11
		nRows = 11
		nCols = 1

		fig =plt.figure(figsize=(3.4,13))
		subplots_to_make = []
		for i in range(1, num_subplots+1):
			subplots_to_make.append((nRows, nCols, i))


		for nrows, ncols, plot_number in subplots_to_make:
			sub = fig.add_subplot(nrows, ncols, plot_number)

			sub.tick_params(which = 'both', direction = 'out', labelsize = 6)

			sub.spines['top'].set_visible(False)
			sub.spines['right'].set_visible(False)
			sub.spines['left'].set_position(('outward', 5))
			if plot_number < (num_subplots + 1 - nCols):
				sub.xaxis.set_ticks_position('none')
				sub.spines['bottom'].set_visible(False)
				sub.set_xticks([])
			else:
				sub.spines['bottom'].set_position(('outward', 5))
		###################

		# Load data from KB
		sim_data = cPickle.load(open(simDataFile, "rb"))
		nAvogadro = sim_data.constants.nAvogadro
		cellDensity = sim_data.constants.cellDensity
		oriC = sim_data.constants.oriCCenter.asNumber()
		terC = sim_data.constants.terCCenter.asNumber()
		genomeLength = len(sim_data.process.replication.genome_sequence)

		recruitmentColNames = sim_data.process.transcription_regulation.recruitmentColNames
		tfs = sorted(set([x.split("__")[-1] for x in recruitmentColNames
			if x.split("__")[-1] != "alpha"]))
		trpRIndex = [i for i, tf in enumerate(tfs) if tf == "CPLX-125"][0]

		tfBoundIds = [target + "__CPLX-125" for target in sim_data.tfToFC["CPLX-125"].keys()]
		synthProbIds = [target + "[c]" for target in sim_data.tfToFC["CPLX-125"].keys()]

		# Initialize for instantaneous growth rate calculations:
		timeMultigen = np.zeros(0)
		cellMassGrowthRateMultigen = np.zeros(0)
		proteinGrowthRateMultigen = np.zeros(0)
		rnaGrowthRateMultigen = np.zeros(0)

		for simDir in allDirs:
			simOutDir = os.path.join(simDir, "simOut")
			# Load time
			initialTime = 0#TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
			time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
			timeStep = TableReader(os.path.join(simOutDir, "Main")).readColumn("timeStepSec")

			# Load mass data
			# Total cell mass is needed to compute concentrations (since we have cell density)
			# Protein mass is needed to compute the mass fraction of the proteome that is trpA
			massReader = TableReader(os.path.join(simOutDir, "Mass"))

			cellMass = massReader.readColumn("cellMass")
			proteinMass = massReader.readColumn("proteinMass")
			rnaMass = massReader.readColumn("rnaMass")
			dnaMass = massReader.readColumn("dnaMass")
			growth_rate = massReader.readColumn("instantaniousGrowthRate")
			cellMass_fg = units.fg * massReader.readColumn("cellMass")
			proteinMass_fg = units.fg * massReader.readColumn("proteinMass")
			massReader.close()

			# Load data from ribosome data listener
			ribosomeDataReader = TableReader(os.path.join(simOutDir, "RibosomeData"))
			nTrpATranslated = ribosomeDataReader.readColumn("numTrpATerminated")
			ribosomeDataReader.close()

			# Load data from bulk molecules
			bulkMoleculesReader = TableReader(os.path.join(simOutDir, "BulkMolecules"))
			bulkMoleculeIds = bulkMoleculesReader.readAttribute("objectNames")

			# Get the concentration of intracellular trp
			trpId = ["TRP[c]"]
			trpIndex = np.array([bulkMoleculeIds.index(x) for x in trpId])
			trpCounts = bulkMoleculesReader.readColumn("counts")[:, trpIndex].reshape(-1)
			trpMols = 1. / nAvogadro * trpCounts
			volume = cellMass_fg / cellDensity
			trpConcentration = trpMols * 1. / volume

			# Get the promoter-bound status for all regulated genes
			tfBoundIndex = np.array([bulkMoleculeIds.index(x) for x in tfBoundIds])
			tfBoundCounts = bulkMoleculesReader.readColumn("counts")[:, tfBoundIndex]

			# Get the amount of monomeric trpA
			trpAProteinId = ["TRYPSYN-APROTEIN[c]"]
			trpAProteinIndex = np.array([bulkMoleculeIds.index(x) for x in trpAProteinId])
			trpAProteinCounts = bulkMoleculesReader.readColumn("counts")[:, trpAProteinIndex].reshape(-1)

			# Get the amount of complexed trpA
			trpABComplexId = ["TRYPSYN[c]"]
			trpABComplexIndex = np.array([bulkMoleculeIds.index(x) for x in trpABComplexId])
			trpABComplexCounts = bulkMoleculesReader.readColumn("counts")[:, trpABComplexIndex].reshape(-1)

			# Get the amount of trpA mRNA
			trpARnaId = ["EG11024_RNA[c]"]
			trpARnaIndex = np.array([bulkMoleculeIds.index(x) for x in trpARnaId])
			trpARnaCounts = bulkMoleculesReader.readColumn("counts")[:, trpARnaIndex].reshape(-1)

			#Get mass per Origin
			mass_per_oriC = TableReader(os.path.join(simOutDir, "ReplicationData")).readColumn("criticalMassPerOriC")

			#RNA over protein:
			#Get active ribosome counts
			uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))
			ribosomeIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRibosome")
			ribosomeCounts = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, ribosomeIndex]
			uniqueMoleculeCounts.close()

			#Find the sequence index and length (to find the fork position later):
			sequenceIdx = TableReader(os.path.join(simOutDir, "ReplicationData")).readColumn("sequenceIdx")
			sequenceLength = TableReader(os.path.join(simOutDir, "ReplicationData")).readColumn("sequenceLength")
			sequenceLength[sequenceLength == -1] = np.nan

			bulkMoleculesReader.close()

			# Computations:

			# Compute total counts and concentration of trpA in monomeric and complexed form
			# (we know the stoichiometry)
			trpAProteinTotalCounts = trpAProteinCounts + 2 * trpABComplexCounts

			# Compute moving averages
			width = 100

			tfBoundCountsMA = np.array([np.convolve(tfBoundCounts[:,i],
				np.ones(width) / width, mode = "same")
					for i in range(tfBoundCounts.shape[1])]).T

			# Find the index of initialization:
			idxInit = np.where(mass_per_oriC >= 1)[0]

			# Calculate the growth rate:
			growth_rate = (1 / units.s) * growth_rate
			growth_rate = growth_rate.asNumber(1 / units.min)

			# Calculate Ribosome Concentration:
			ribosomeConcentration = ((1 / sim_data.constants.nAvogadro) * ribosomeCounts) / ((1.0 / sim_data.constants.cellDensity) * (cellMass_fg))
			ribosomeConcentration = ribosomeConcentration.asNumber(units.umol / units.L)

			# Fork Position:
			reverseIdx = 1
			reverseCompIdx = 3
			reverseSequences = np.logical_or(sequenceIdx == reverseIdx,
				sequenceIdx == reverseCompIdx)
			sequenceLength[reverseSequences] = -1 * sequenceLength[reverseSequences]

			# Down sample dna polymerase position, every position is only plotted once here
			unique, index, value = np.unique(sequenceLength, return_index=True,
				return_inverse=True)
			m = np.zeros_like(value, dtype=bool)
			m[index] = True
			m = m.reshape(sequenceLength.shape)
			sequenceLength[~m] = np.nan

			dnap_times = []
			dnap_positions = []
			for t, pos in zip(time, sequenceLength)[::10]:
				for p in pos[np.isfinite(pos)]:
					dnap_times += [t]
					dnap_positions += [p]
			dnap_times = np.array(dnap_times)
			dnap_positions = np.array(dnap_positions)

			# Relative Rate of dNTP polymerization
			relative_rate_dNTP_poly = (sequenceIdx != -1).sum(axis = 1) / 4

			# Calculate the average instananeous growth rate

			cellmassGrowthRate = np.diff(cellMass) / cellMass[1:] / timeStep[:-1]
			proteinGrowthRate = np.diff(proteinMass) / proteinMass[1:] / timeStep[:-1]
			rnaGrowthRate = np.diff(rnaMass) / rnaMass[1:] / timeStep[:-1]
			dnaGrowthRate = np.diff(dnaMass) / dnaMass[1:] / timeStep[:-1]

			timeMultigen = np.hstack((timeMultigen, time[:-1]))
			cellMassGrowthRateMultigen = np.hstack((cellMassGrowthRateMultigen, cellmassGrowthRate))
			proteinGrowthRateMultigen = np.hstack((proteinGrowthRateMultigen, proteinGrowthRate))
			rnaGrowthRateMultigen = np.hstack((rnaGrowthRateMultigen, rnaGrowthRate))

			'''
			Plots:
			Note: Some of the Y axis limits are hard coded. 
			This was to make the plot have more rounded numbers for the paper. 
			These limits may change simulation to simulation though, 
			so make sure to change to be more dynamic 
			if using on a different data set.
			'''

			# Plot parameters:
			plot_line_color = '#0d71b9'
			plot_marker_color = '#ed2224'
			plot_font_size = 6

			##############################################################
			# Plot Cell Mass
			ax1 = plt.subplot(nRows, nCols, 1)
			ax1.plot(time, cellMass, color = plot_line_color)
			ax1.plot(time[idxInit], cellMass[idxInit],  markersize=4,
				linewidth=0, marker="o", color = plot_marker_color,
				markeredgewidth=0)
			plt.ylabel("Cell Mass\n(fg)", fontsize = plot_font_size)
			y_min_1, y_max_1 = 1000, 4600
			ax1.set_ylim([y_min_1, y_max_1])
			ax1.set_yticks([y_min_1, y_max_1])
			ax1.set_xlim([0, time.max()])
			##############################################################
			# Plot Instantaneous Growth Rate
			ax2 = plt.subplot(nRows, nCols, 2)
			ax2.plot(time, growth_rate, color = plot_line_color)
			plt.ylabel("Instantaneouse growth Rate", fontsize = plot_font_size)
			ax2.set_ylabel(r"$\mu$ $(\frac{gDCW}{gDCW \cdot \, min})$")
			y_min_2, y_max_2 = 0, 0.032
			ax2.set_ylim([y_min_2, y_max_2])
			ax2.set_yticks([y_min_2, y_max_2])
			ax2.set_xlim([0, time.max()])
			##############################################################
			# Plot Active Ribosome Concentration
			ax3 = plt.subplot(nRows, nCols, 3)
			ax3.plot(time, ribosomeConcentration, color = plot_line_color)
			ax3.set_ylim([15., 25.])
			plt.ylabel("Active Ribosome\n(umol/L)", fontsize = plot_font_size)
			y_min_3, y_max_3 = ax3.get_ylim()
			ax3.set_yticks([y_min_3, y_max_3])
			ax3.set_yticklabels(["%0.0f" % y_min_3, "%0.0f" % y_max_3])
			ax3.set_xlim([0, time.max()])
			##############################################################
			# Plot DNA Polymerase Position
			ax4 = plt.subplot(nRows, nCols, 4)
			ax4.plot(dnap_times, dnap_positions, marker='o', markersize=1,
				linewidth=0, color = plot_line_color)
			plt.ylabel("DNA polymerase\nposition", fontsize = plot_font_size)
			ax4.set_yticks([-1 * genomeLength / 2, 0, genomeLength / 2])
			ax4.set_yticklabels(['-terC', 'oriC', '+terC'])
			ax4.set_xlim([0, time.max()])
			##############################################################
			#Plot Relative Rate of dNTP Polymerization
			ax5 = plt.subplot(nRows, nCols, 5)
			ax5.plot(time, relative_rate_dNTP_poly, color = plot_line_color)
			plt.ylabel("Relative rate of\ndNTP polymerization",
				fontsize = plot_font_size)
			y_min_5, y_max_5 = 0, 6.
			ax5.set_ylim([y_min_5, y_max_5])
			ax5.set_yticks([y_min_5, y_max_5])
			ax5.set_xlim([0, time.max()])
			##############################################################
			#Plot Promoter Occupancy
			ax6 = plt.subplot(nRows, nCols, 6)
			ax6.plot(time, tfBoundCountsMA, color = plot_line_color)
			#comment out color in order to see colors per generation
			plt.ylabel("TrpR Bound To Promoters\n(Moving Average)",
				fontsize = plot_font_size)
			y_min_6, y_max_6 = 0, 1
			ax6.set_ylim([y_min_6, y_max_6])
			ax6.set_yticks([y_min_6, y_max_6])
			ax6.set_xlim([0, time.max()])
			##############################################################
			#Plot TrpA mRNA Counts
			ax7 = plt.subplot(nRows, nCols, 7)
			ax7.plot(time, trpARnaCounts, color = plot_line_color)
			plt.ylabel("TrpA mRNA\nCounts", fontsize = plot_font_size)
			y_min_7, y_max_7 = ax7.get_ylim()
			ax7.set_ylim([0, y_max_7])
			ax7.set_yticks([0, y_max_7])
			ax7.set_yticklabels([0, "%0.0f" % y_max_7])
			ax7.set_xlim([0, time.max()])
			##############################################################
			# Plot Translation Events
			ax8 = plt.subplot(nRows, nCols, 8)
			ax8.plot(time, nTrpATranslated, color = '#0d71b9')
			plt.ylabel("Number of TrpA translation events",
				fontsize = plot_font_size)
			y_min_8, y_max_8 = ax8.get_ylim()
			ax8.set_yticks([y_min_8, y_max_8])
			ax8.set_yticklabels(["%0.0f" % y_min_8, "%0.0f" % y_max_8])
			ax8.set_xlim([0, time.max()])
			##############################################################
			# Plot Total TrpA Counts
			ax9 = plt.subplot(nRows, nCols, 9)
			ax9.plot(time, trpAProteinTotalCounts, color = plot_line_color)
			plt.ylabel("TrpA Counts", fontsize = plot_font_size)
			y_min_9, y_max_9 = 800, 6000.
			ax9.set_ylim([y_min_9, y_max_9])
			ax9.set_yticks([y_min_9,y_max_9])
			ax9.set_xlim([0, time.max()])
			##############################################################
			# Plot Internal TRP concentraion
			ax10 = plt.subplot(nRows, nCols, 10)
			ax10.plot(time, trpConcentration.asNumber(units.umol / units.L),
				color = plot_line_color)
			plt.ylabel("Internal TRP Conc.\n(uM)", fontsize = plot_font_size)
			y_min_10, y_max_10 = 0, 400
			ax10.set_ylim([y_min_10, y_max_10])
			ax10.set_yticks([y_min_10, y_max_10])
			ax10.set_yticklabels([y_min_10, y_max_10])
			ax10.set_xlim([0, time.max()])
			##############################################################

		# Get the moving averages of the instantaneous growth rates.

		width_inst = 200
		cellMassGrowthRateMultigen = np.convolve(cellMassGrowthRateMultigen, np.ones(width_inst) / width_inst, mode = "same")
		proteinGrowthRateMultigen = np.convolve(proteinGrowthRateMultigen, np.ones(width_inst) / width_inst, mode = "same")
		rnaGrowthRateMultigen = np.convolve(rnaGrowthRateMultigen, np.ones(width_inst) / width_inst, mode = "same")
		##############################################################
		# Plot average instantanous growth rate
		colors = ["#43aa98", "#0071bb", "#bf673c"]
		ax11 = plt.subplot(nRows, nCols, 11)
		ax11.plot(timeMultigen[:-width_inst], cellMassGrowthRateMultigen[:-width_inst] * 60.,color = colors[0])
		ax11.plot(timeMultigen[:-width_inst], proteinGrowthRateMultigen[:-width_inst] * 60., color = colors[1])
		ax11.plot(timeMultigen[:-width_inst], rnaGrowthRateMultigen[:-width_inst] * 60., color = colors[2])

		plt.ylabel("Average Instantaneous Growth Rate.\n(min -1)", fontsize = plot_font_size)
		y_min_11, y_max_11 = 0.014, 0.029
		ax11.set_ylim([y_min_11, y_max_11])
		ax11.set_yticks([y_min_11, y_max_11])
		ax11.set_yticklabels([y_min_11, y_max_11])
		ax11.set_xlim([0, time.max()])
		##############################################################

		for nrows, ncols, plot_number in subplots_to_make:
			if plot_number > (num_subplots - nCols):
				sub.set_xticks([0, time.max()])
				sub.set_xticklabels([0., np.round(time.max() / 60.,
					decimals = 0)])
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
