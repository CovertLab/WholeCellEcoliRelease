"""
Plots time-step effects on polymerization

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/18/2015
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
		rnaIds = sim_data.process.transcription.rnaData['id']

		elongationRate = float(sim_data.constants.ribosomeElongationRateMax.asNumber(units.aa / units.s))
		rnaLengths = sim_data.process.transcription.rnaData['length'].asNumber()

		rRnaIds = []
		rRnaIds.append(sim_data.moleculeGroups.s30_16sRRNA[0]) # Only using the rrnA operon subunits! This is a bit of a hack...
		rRnaIds.append(sim_data.moleculeGroups.s50_23sRRNA[0])
		rRnaIds.append(sim_data.moleculeGroups.s50_5sRRNA[0])

		rRnaIndexes = np.array([np.where(rnaIds == comp)[0][0] for comp in rRnaIds])
		rRnaLengths = rnaLengths[rRnaIndexes]

		polymerizedNtpsInActiveRibosome = rRnaLengths.sum()

		# Load time
		initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
		timeStep = units.s * TableReader(os.path.join(simOutDir, "Main")).readColumn("timeStepSec")

		# Load ribosome data
		ribosomeDataFile = TableReader(os.path.join(simOutDir, "RibosomeData"))
		ribosomesTerminated = ribosomeDataFile.readColumn("didTerminate")
		ribosomesInitialized = ribosomeDataFile.readColumn("didInitialize")
		ribosomeTerminationLoss = ribosomeDataFile.readColumn("terminationLoss")
		actualRibosomeElongations = ribosomeDataFile.readColumn("actualElongations")
		ribosomeDataFile.close()

		# Load RNAP data
		rnapDataFile = TableReader(os.path.join(simOutDir, "RnapData"))
		rnapsTerminated = rnapDataFile.readColumn("didTerminate")
		rnapsInitialized = rnapDataFile.readColumn("didInitialize")
		rnapTerminationLoss = rnapDataFile.readColumn("terminationLoss")
		actualRnapElongations = rnapDataFile.readColumn("actualElongations")
		rnapDataFile.close()

		# Load count data for s30 proteins, rRNA, and final 30S complex
		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		bulkMoleculeCounts = bulkMolecules.readColumn("counts")

		# Get indexes
		moleculeIds = bulkMolecules.readAttribute("objectNames")
		ribosomeSubunitIndexes = np.array([moleculeIds.index(comp) for comp in ribosomeSubunitIds], np.int)
		rnaIndexes = np.array([moleculeIds.index(comp) for comp in rnaIds])

		# Load data
		ribosomeSubunitCounts = bulkMoleculeCounts[:, ribosomeSubunitIndexes]
		rnaCounts = bulkMoleculeCounts[:, rnaIndexes]

		bulkMolecules.close()

		uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))

		ribosomeIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRibosome")
		activeRibosome = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, ribosomeIndex]

		uniqueMoleculeCounts.close()

		rnaDegradationListenerFile = TableReader(os.path.join(simOutDir, "RnaDegradationListener"))
		countRnaDegraded = rnaDegradationListenerFile.readColumn('countRnaDegraded')
		rnaDegradationListenerFile.close()

		# Calculate statistics
		totalRibosome = (activeRibosome + np.min(ribosomeSubunitCounts))
		totalRibosomeCapacity = totalRibosome * elongationRate

		freeSubunitMass = (ribosomeSubunitMasses * ribosomeSubunitCounts / nAvogadro).asNumber(units.fg).sum(axis = 1)
		activeRibosomeMass = (mass70s * activeRibosome / nAvogadro).asNumber(units.fg)
		totalRibosomeMass = freeSubunitMass + activeRibosomeMass
		massFractionActive = activeRibosomeMass / totalRibosomeMass

		bulkPolymerizedNtpCounts = np.dot(rnaCounts, rnaLengths)
		uniquePolymerizedNtpCounts = activeRibosome * polymerizedNtpsInActiveRibosome
		totalPolymerizedNtpCounts = bulkPolymerizedNtpCounts + uniquePolymerizedNtpCounts
		polymerizedNtpsDegraded = np.dot(countRnaDegraded, rnaLengths)

		growthRate = np.log(2) / sim_data.doubling_time
		expectedGrowth = totalPolymerizedNtpCounts * (np.exp(units.convertNoUnitToNumber(growthRate * timeStep)) - 1)
		actualGrowth = actualRnapElongations - polymerizedNtpsDegraded

		plt.figure(figsize = (8.5, 11))

		ribosomeInit_axis = plt.subplot(9,1,1)
		ribosomeInit_axis.plot(time / 60., ribosomesInitialized, label="Number of ribosomes initialized", linewidth=2)
		ribosomeInit_axis.set_ylabel("Ribosomes\ninitialized")

		ribosomeTerm_axis = plt.subplot(9,1,2)
		ribosomeTerm_axis.plot(time / 60., ribosomesTerminated, label="Number of ribosomes terminated", linewidth=2)
		ribosomeTerm_axis.set_ylabel("Ribosomes\nterminated")

		ribosomeTermLoss_axis = plt.subplot(9,1,3)
		ribosomeTermLoss_axis.plot(time / 60., ribosomeTerminationLoss, label="Lost capacity due to termination", linewidth=2)
		ribosomeTermLoss_axis.set_ylabel("Ribosome\ntermination\nloss")

		ribosomePercentLoss_axis = plt.subplot(9,1,4)
		ribosomePercentLoss_axis.plot(time / 60., ribosomeTerminationLoss / actualRibosomeElongations.astype(np.float) * 100, label="Lost capacity due to termination", linewidth=2)
		ribosomePercentLoss_axis.set_ylabel("Ribosome percent\nloss")

		rnapInit_axis = plt.subplot(9,1,5)
		rnapInit_axis.plot(time / 60., rnapsInitialized, label="Number of rnap initialized", linewidth=2)
		rnapInit_axis.set_ylabel("Rnap\ninitalized")

		rnapTerm_axis = plt.subplot(9,1,6)
		rnapTerm_axis.plot(time / 60., rnapsTerminated, label="Number of rnap terminated", linewidth=2)
		rnapTerm_axis.set_ylabel("Rnap\nterminated")

		rnapTermLoss_axis = plt.subplot(9,1,7)
		rnapTermLoss_axis.plot(time / 60., rnapTerminationLoss, label="Lost capacity due to termination", linewidth=2)
		rnapTermLoss_axis.set_ylabel("RNAP\ntermination\nloss")

		rnapActualVsExpectedGrowth_axis = plt.subplot(9,1,8)
		rnapActualVsExpectedGrowth_axis.plot(time / 60., actualGrowth / expectedGrowth.astype(np.float), label="Lost capacity due to termination", linewidth=2)
		rnapActualVsExpectedGrowth_axis.set_ylim([0, 2])
		rnapActualVsExpectedGrowth_axis.set_ylabel("Actual growth/\nexpected")

		# Save
		plt.subplots_adjust(hspace = 0.5, wspace = 0.5)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
