"""
Plots ribosome capacity

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 11/20/2014
"""

from __future__ import absolute_import

import os

import numpy as np
from matplotlib import pyplot as plt
import cPickle

from wholecell.io.tablereader import TableReader
from wholecell.utils import units
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot

FONT = {
		'size'	:	8
		}


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Load data from KB
		sim_data = cPickle.load(open(simDataFile, "rb"))
		nAvogadro = sim_data.constants.nAvogadro
		ribosomeSubunitIds = []
		ribosomeSubunitIds.append(sim_data.moleculeIds.s50_fullComplex)
		ribosomeSubunitIds.append(sim_data.moleculeIds.s30_fullComplex)
		ribosomeSubunitIds.extend(sim_data.moleculeGroups.s50_proteinComplexes)
		ribosomeSubunitIds.extend(sim_data.process.complexation.getMonomers(sim_data.moleculeIds.s50_fullComplex)['subunitIds'])
		ribosomeSubunitIds.extend(sim_data.process.complexation.getMonomers(sim_data.moleculeIds.s30_fullComplex)['subunitIds'])
		ribosomeSubunitMasses = sim_data.getter.getMass(ribosomeSubunitIds)
		mass70s = (sim_data.getter.getMass([sim_data.moleculeIds.s50_fullComplex]) + sim_data.getter.getMass([sim_data.moleculeIds.s30_fullComplex]))[0]

		timeStep = TableReader(os.path.join(simOutDir, "Main")).readColumn("timeStepSec")
		elongationRate = float(sim_data.growthRateParameters.ribosomeElongationRate.asNumber(units.aa / units.s)) * timeStep

		# Load ribosome data
		ribosomeDataFile = TableReader(os.path.join(simOutDir, "RibosomeData"))

		actualElongations = ribosomeDataFile.readColumn("actualElongations")
		expectedElongations_recorded = ribosomeDataFile.readColumn("expectedElongations")
		initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

		ribosomeDataFile.close()

		# Load count data for s30 proteins, rRNA, and final 30S complex
		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))

		# Get indexes
		moleculeIds = bulkMolecules.readAttribute("objectNames")
		ribosomeSubunitIndexes = np.array([moleculeIds.index(comp) for comp in ribosomeSubunitIds], np.int)

		# Load data
		ribosomeSubunitCounts = bulkMolecules.readColumn("counts")[:, ribosomeSubunitIndexes]

		bulkMolecules.close()

		uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))

		ribosomeIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRibosome")
		activeRibosome = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, ribosomeIndex]

		uniqueMoleculeCounts.close()

		# Calculate statistics
		totalRibosome = (activeRibosome + ribosomeSubunitCounts.min(axis=1))
		totalRibosomeCapacity = totalRibosome * elongationRate

		freeSubunitMass = (ribosomeSubunitMasses * ribosomeSubunitCounts / nAvogadro).asNumber(units.fg).sum(axis = 1)
		activeRibosomeMass = (mass70s * activeRibosome / nAvogadro).asNumber(units.fg)
		totalRibosomeMass = freeSubunitMass + activeRibosomeMass
		massFractionActive = activeRibosomeMass / totalRibosomeMass

		plt.figure(figsize = (8.5, 11))
		plt.rc('font', **FONT)

		ribosomeCapacity_axis = plt.subplot(6,1,1)
		ribosomeCapacity_axis.plot(time / 60., totalRibosomeCapacity, label="Theoretical total ribosome capacity", linewidth=2, color='b')
		ribosomeCapacity_axis.plot(time / 60., actualElongations, label="Actual elongations", linewidth=2, color='r')
		ribosomeCapacity_axis.set_ylabel("Amino acids\npolymerized")
		ribosomeCapacity_axis.legend(ncol=2)

		activeRibosomeCapacity_axis = plt.subplot(6,1,2)
		activeRibosomeCapacity_axis.plot(time / 60., activeRibosome * elongationRate, label="Theoretical active ribosome capacity", linewidth=2, color='b')
		activeRibosomeCapacity_axis.plot(time / 60., actualElongations, label="Actual elongations", linewidth=2, color='r')
		activeRibosomeCapacity_axis.set_ylabel("Amino acids\npolymerized")
		activeRibosomeCapacity_axis.legend(ncol=2)

		inactiveRibosomeCapacity_axis = plt.subplot(6,1,3)
		inactiveRibosomeCapacity_axis.plot(time / 60., ribosomeSubunitCounts.min(axis=1) * elongationRate, label="Theoretical inactive ribosome capacity", linewidth=2, color='b')
		inactiveRibosomeCapacity_axis.set_ylabel("Amino acids\npolymerized")
		inactiveRibosomeCapacity_axis.legend(ncol=2)

		fractionalCapacity_axis = plt.subplot(6,1,4)
		fractionalCapacity_axis.plot(time / 60., actualElongations / totalRibosomeCapacity, label="Fraction of ribosome capacity used", linewidth=2, color='k')
		fractionalCapacity_axis.set_ylabel("Fraction of\nribosome capacity\nused")
		fractionalCapacity_axis.set_yticks(np.arange(0., 1.05, 0.1))
		#fractionalCapacity_axis.get_yaxis().grid(b=True, which='major', color='b', linestyle='--')
		fractionalCapacity_axis.grid(b=True, which='major', color='b', linestyle='--')

		effectiveElongationRate_axis = plt.subplot(6,1,5)
		effectiveElongationRate_axis.plot(time / 60., actualElongations / activeRibosome, label="Effective elongation rate", linewidth=2, color='k')
		effectiveElongationRate_axis.set_ylabel("Effective\nelongation rate\n(aa/s/ribosome)")

		fractionActive_axis = plt.subplot(6,1,6)
		fractionActive_axis.plot(time / 60., massFractionActive, label="Mass fraction active", linewidth=2, color='k')
		fractionActive_axis.set_ylabel("Mass fraction of\nactive ribosomes")
		fractionActive_axis.set_yticks(np.arange(0., 1.1, 0.1))
		fractionActive_axis.grid(b=True, which='major', color='b', linestyle='--')

		# Save
		plt.subplots_adjust(hspace = 0.5, wspace = 0.6, top = 0.95, bottom = 0.05)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
