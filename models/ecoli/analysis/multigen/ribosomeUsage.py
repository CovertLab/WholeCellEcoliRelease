"""
Plot usage statistics of ribosomes

@author: Gwanggyu Sun
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 10/18/2017
"""

from __future__ import absolute_import, division

import os
import cPickle

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.utils import units
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(seedOutDir):
			raise Exception, "seedOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)

		# Get first cell from each generation
		firstCellLineage = []

		# For all generation indexes subject to analysis, get first cell
		for gen_idx in range(ap.n_generation):
			firstCellLineage.append(ap.get_cells(generation = [gen_idx])[0])

		# Get sim data from cPickle file
		sim_data = cPickle.load(open(simDataFile, "rb"))

		# Create new figure and set size
		fig = plt.figure()
		fig.set_size_inches(15,12)

		# Divide figure into subplot grids
		gs = gridspec.GridSpec(7, 2)

		ax1  = plt.subplot(gs[0,0])
		ax2  = plt.subplot(gs[1,0])
		ax3  = plt.subplot(gs[2,0])
		ax4  = plt.subplot(gs[3,0])
		ax5  = plt.subplot(gs[4,0])
		ax6  = plt.subplot(gs[5,0])
		ax7  = plt.subplot(gs[6,0])

		ax8  = plt.subplot(gs[0,1])
		ax9  = plt.subplot(gs[1,1])
		ax10 = plt.subplot(gs[2,1])
		ax11 = plt.subplot(gs[3,1])
		ax12 = plt.subplot(gs[4,1])
		ax13 = plt.subplot(gs[5,1])
		ax14 = plt.subplot(gs[6,1])


		# Go through first cells in each generation
		for gen, simDir in enumerate(firstCellLineage):

			simOutDir = os.path.join(simDir, "simOut")

			## Mass growth rate ##
			time, growthRate = getMassData(simDir, ["instantaniousGrowthRate"])
			timeStep = units.s * TableReader(os.path.join(simOutDir, "Main")).readColumn("timeStepSec")
			time = units.s * time

			## Ribosome counts and statistics ##

			# Get ids for 30S and 50S subunits
			complexIds30S = [sim_data.moleculeIds.s30_fullComplex]
			complexIds50S = [sim_data.moleculeIds.s50_fullComplex]

			# Get molecular weights for 30S and 50S subunits, and add these two for 70S
			nAvogadro = sim_data.constants.nAvogadro
			mw30S = sim_data.getter.getMass(complexIds30S)
			mw50S = sim_data.getter.getMass(complexIds50S)
			mw70S = mw30S + mw50S

			# Get indexes for 30S and 50S subunits based on ids
			bulkMoleculesDataFile = TableReader(os.path.join(simOutDir, "BulkMolecules"))
			moleculeIds = bulkMoleculesDataFile.readAttribute("objectNames")
			bulkMoleculeCounts = bulkMoleculesDataFile.readColumn("counts")

			complexIndexes30S = np.array([moleculeIds.index(comp) for comp in complexIds30S], np.int)
			complexIndexes50S = np.array([moleculeIds.index(comp) for comp in complexIds50S], np.int)

			# Get counts of 30S and 50S mRNA, rProteins, rRNA, and full complex counts
			complexCounts30S = bulkMoleculeCounts[:, complexIndexes30S]
			complexCounts50S = bulkMoleculeCounts[:, complexIndexes50S]

			bulkMoleculesDataFile.close()

			# Get active ribosome counts
			uniqueMoleculeCountsDataFile = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))

			ribosomeIndex = uniqueMoleculeCountsDataFile.readAttribute("uniqueMoleculeIds").index("activeRibosome")
			activeRibosome = uniqueMoleculeCountsDataFile.readColumn("uniqueMoleculeCounts")[:, ribosomeIndex]

			uniqueMoleculeCountsDataFile.close()

			# Get ribosome data
			ribosomeDataFile = TableReader(os.path.join(simOutDir, "RibosomeData"))

			didInitialize = ribosomeDataFile.readColumn("didInitialize")
			actualElongations = ribosomeDataFile.readColumn("actualElongations")
			didTerminate = ribosomeDataFile.readColumn("didTerminate")
			effectiveElongationRate = ribosomeDataFile.readColumn("effectiveElongationRate")

			ribosomeDataFile.close()

			# Get mass data
			massDataFile = TableReader(os.path.join(simOutDir, "Mass"))

			cellMass = massDataFile.readColumn("cellMass")

			massDataFile.close()

			# Calculate cell volume
			cellVolume = (1.0 / sim_data.constants.cellDensity) * (units.fg * cellMass)

			# Calculate molecule counts and molar fraction of active ribosomes
			counts30S = complexCounts30S
			counts50S = complexCounts50S
			activeRibosomeCounts = activeRibosome
			totalRibosomeCounts = activeRibosomeCounts + np.hstack((counts30S, counts50S)).min(axis=1)
			molarFractionActive = activeRibosomeCounts.astype(np.float) / totalRibosomeCounts

			totalRibosomeConcentration = ((1 / nAvogadro) * totalRibosomeCounts) / cellVolume
			activeRibosomeConcentration = ((1 / nAvogadro) * activeRibosomeCounts) / cellVolume

			# Calculate molecule masses and mass fraction of active ribosomes
			mass30S = ((1 / nAvogadro) * counts30S) * mw30S
			mass50S = ((1 / nAvogadro) * counts50S) * mw50S
			activeRibosomeMass = ((1 / nAvogadro) * np.transpose([activeRibosomeCounts])) * mw70S

			totalRibosomeMass   = activeRibosomeMass + mass30S + mass50S
			massFractionActive  = activeRibosomeMass / totalRibosomeMass

			# ax1: Plot timestep
			ax1.plot(time.asNumber(units.min), timeStep.asNumber(units.s), linestyle='-')
			if gen == ap.n_generation - 1:
				ax1.set_xlim([-5, max(time.asNumber(units.min))])
				ax1.set_ylim([0, 1])
			ax1.set_ylabel("Length of\ntime step (s)")

			# ax2: Plot cell volume
			ax2.plot(time.asNumber(units.min), cellVolume.asNumber(units.L), linestyle='-')
			if gen == ap.n_generation - 1:
				ax2.set_xlim([-5, max(time.asNumber(units.min))])
				ax2.set_ylim([0, 3e-15])
			ax2.set_ylabel("Cell volume\n(L)")

			# ax3: Plot total ribosome counts
			ax3.plot(time.asNumber(units.min), totalRibosomeCounts, linestyle='-')
			if gen == ap.n_generation - 1:
				ax3.set_xlim([-5, max(time.asNumber(units.min))])
				ax3.set_ylim([10000, 35000])
			ax3.set_ylabel("Total ribosome\ncount")

			# ax4: Plot total ribosome concentrations
			ax4.plot(time.asNumber(units.min), totalRibosomeConcentration.asNumber(units.mmol / units.L), linestyle='-')
			if gen == ap.n_generation - 1:
				ax4.set_xlim([-5, max(time.asNumber(units.min))])
				ax4.set_ylim([0.019, 0.023])
			ax4.set_ylabel("[Total ribosome]\n(mM)")

			# ax5: Plot active ribosome counts
			if gen == 0:
				ax5.plot(time[1:].asNumber(units.min), activeRibosomeCounts[1:], linestyle='-')
			else:
				ax5.plot(time.asNumber(units.min), activeRibosomeCounts, linestyle='-')
			if gen == ap.n_generation - 1:
				ax5.set_xlim([-5, max(time.asNumber(units.min))])
				ax5.set_ylim([10000, 30000])
			ax5.set_ylabel("Active ribosome\ncount")

			# ax6: Plot active ribosome concentrations
			if gen == 0:
				ax6.plot(time[1:].asNumber(units.min), activeRibosomeConcentration[1:].asNumber(units.mmol / units.L), linestyle='-')
			else:
				ax6.plot(time.asNumber(units.min), activeRibosomeConcentration.asNumber(units.mmol / units.L), linestyle='-')
			if gen == ap.n_generation - 1:
				ax6.set_xlim([-5, max(time.asNumber(units.min))])
				ax6.set_ylim([0.0, 0.023])
			ax6.set_ylabel("[Active ribosome]\n(mM)")

			# ax7: Plot molar fraction of active ribosomes
			if gen == 0:
				ax7.plot(time[1:].asNumber(units.min), molarFractionActive[1:], linestyle='-')
			else:
				ax7.plot(time.asNumber(units.min), molarFractionActive, linestyle='-')
			if gen == ap.n_generation - 1:
				ax7.set_xlim([-5, max(time.asNumber(units.min))])
				ax7.set_ylim([0.7, 1])
			ax7.set_ylabel("Molar fraction\nactive ribosomes")

			# ax8: Plot mass fraction of active ribosomes
			if gen == 0:
				ax8.plot(time[1:].asNumber(units.min), massFractionActive[1:], linestyle='-')
			else:
				ax8.plot(time.asNumber(units.min), massFractionActive, linestyle='-')
			if gen == ap.n_generation - 1:
				ax8.set_xlim([-5, max(time.asNumber(units.min))])
				ax8.set_ylim([0.7, 1])
			ax8.set_ylabel("Mass fraction\nactive ribosomes")

			# ax9: Plot number of activations
			ax9.plot(time.asNumber(units.min), didInitialize, linestyle='-')
			if gen == ap.n_generation - 1:
				ax9.set_xlim([-5, max(time.asNumber(units.min))])
				ax9.set_ylim([0, 2000])
			ax9.set_ylabel("Activations\nper timestep")

			# ax10: Plot number of deactivations (terminated translations)
			ax10.plot(time.asNumber(units.min), didTerminate, linestyle='-')
			if gen == ap.n_generation - 1:
				ax10.set_xlim([-5, max(time.asNumber(units.min))])
				ax10.set_ylim([0, 2000])
			ax10.set_ylabel("Deactivations\nper timestep")

			# ax11: Plot number of activations per time * volume
			ax11.plot(time.asNumber(units.min), didInitialize / (timeStep.asNumber(units.s) * cellVolume.asNumber(units.L)), linestyle='-')
			if gen == ap.n_generation - 1:
				ax11.set_xlim([-5, max(time.asNumber(units.min))])
				ax11.set_ylim([0, 1.2e18])
			ax11.set_ylabel("Activations\nper time*volume")

			# ax12: Plot number of deactivations per time * volume
			ax12.plot(time.asNumber(units.min), didTerminate / (timeStep.asNumber(units.s) * cellVolume.asNumber(units.L)), linestyle='-')
			if gen == ap.n_generation - 1:
				ax12.set_xlim([-5, max(time.asNumber(units.min))])
				ax12.set_ylim([0, 1.2e18])
			ax12.set_ylabel("Deactivations\nper time*volume")

			# ax13: Plot number of amino acids translated in each timestep
			ax13.plot(time.asNumber(units.min), actualElongations, linestyle='-')
			if gen == ap.n_generation - 1:
				ax13.set_xlim([-5, max(time.asNumber(units.min))])
				# ax13.set_ylim([0])
			ax13.set_ylabel("AA translated")

			# ax14: Plot effective ribosome elongation rate for each timestep
			ax14.plot(time.asNumber(units.min), effectiveElongationRate, linestyle='-')
			if gen == ap.n_generation - 1:
				ax14.set_xlim([-5, max(time.asNumber(units.min))])
				ax14.set_ylim([10, 22])
			ax14.set_ylabel("Effective\nelongation rate")

		fig.subplots_adjust(hspace = 0.5, wspace = 0.3)

		exportFigure(plt, plotOutDir, plotOutFileName,metadata)
		plt.close("all")


def getMassData(simDir, massNames):
	simOutDir = os.path.join(simDir, "simOut")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
	mass = TableReader(os.path.join(simOutDir, "Mass"))

	massFractionData = np.zeros((len(massNames),time.size))

	for idx, massType in enumerate(massNames):
		massFractionData[idx,:] = mass.readColumn(massNames[idx])

	if len(massNames) == 1:
		massFractionData = massFractionData.reshape(-1)

	return time, massFractionData


if __name__ == "__main__":
	Plot().cli()
