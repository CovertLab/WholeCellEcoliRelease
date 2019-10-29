"""
Plot mRNA half lives (observed vs. actual)

@author: Javier Carrera
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/30/2015
"""

from __future__ import absolute_import
from __future__ import division

import os

import numpy as np
from matplotlib import pyplot as plt
import cPickle

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.analysis.analysis_tools import read_bulk_molecule_counts
from models.ecoli.analysis import singleAnalysisPlot


# TODO: account for complexation

class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Get the names of rnas from the KB

		sim_data = cPickle.load(open(simDataFile))

		isMRna = sim_data.process.transcription.rnaData["isMRna"]
		rnaIds = sim_data.process.transcription.rnaData["id"][isMRna]

		expectedDegradationRate = sim_data.process.transcription.rnaData['degRate'][isMRna].asNumber()

		rnaCountsBulk, = read_bulk_molecule_counts(simOutDir, rnaIds)
		rnaCounts = rnaCountsBulk[1:,:]
		rnaCountsTotal = rnaCounts.sum(axis = 0)

		rnaDegradationListenerFile = TableReader(os.path.join(simOutDir, "RnaDegradationListener"))
		countRnaDegraded = rnaDegradationListenerFile.readColumn('countRnaDegraded')
		rnaDegradationListenerFile.close()
		rnaDegraded = countRnaDegraded[1:,:]
		rnaDegradedTotal = rnaDegraded.sum(axis = 0)[isMRna]

		rnaSynthesizedListenerFile = TableReader(os.path.join(simOutDir, "TranscriptElongationListener"))
		countRnaSynthesized = rnaSynthesizedListenerFile.readColumn('countRnaSynthesized')
		rnaSynthesizedListenerFile.close()
		rnaSynthesized = countRnaSynthesized[1:,:]
		rnaSynthesizedTotal = rnaSynthesized.sum(axis = 0)[isMRna]

		rnaDegradationRate1 = []
		rnaDegradationRate2 = []
		rnaDegradationRate3 = []
		rnaDegradationRate4 = []
		rnaDegradationRate5 = []
		expectedDegradationRateSubset = []
		expectedDegradationRateSubset4 = []
		expectedDegradationRateSubset5 = []

		for i in range(0, len(rnaCountsTotal)): # loop for genes
			if rnaCountsTotal[i] != 0:
				rnaDegradationRate1.append(rnaDegradedTotal[i] / rnaCountsTotal[i]) # Sum_tau(kd*r) / Sum_tau(r)

				rnaDegradationRate2.append(rnaSynthesizedTotal[i] / rnaCountsTotal[i]) # Sum_tau(sim_data) / Sum_tau(r)

				rnaDegradationRate3.append( (rnaSynthesizedTotal[i] - rnaCounts[0,i]) / rnaCountsTotal[i]) # (Sum_tau(sim_data) - r) / Sum_tau(r)

				rnaDegradationRate4_t = []
				rnaDegradationRate5_t = []
				for j in range(0, len(rnaDegraded[:,i]) - 1): # loop for timepoints
					if rnaCounts[j,i] != 0:
						if rnaDegraded[j,i] > 0:
							rnaDegradationRate4_t.append(rnaDegraded[j,i] / rnaCounts[j,i]) # Average_t(kd*r / r)
						if (rnaSynthesized[j,i] - (rnaCounts[j+1,i] - rnaCounts[j,i])) > 0:
							rnaDegradationRate5_t.append((rnaSynthesized[j,i] - (rnaCounts[j+1,i] - rnaCounts[j,i])) / rnaCounts[j,i]) # Average_t(kd*r / r)

				if len(rnaDegradationRate4_t) > 0:
					rnaDegradationRate4.append(np.mean(rnaDegradationRate4_t))
					expectedDegradationRateSubset4.append(expectedDegradationRate[i])
				if len(rnaDegradationRate4_t) <= 0:
					rnaDegradationRate4.append(-1)
					expectedDegradationRateSubset4.append(-1)


				if len(rnaDegradationRate5_t) > 0:
					rnaDegradationRate5.append(np.mean(rnaDegradationRate5_t))
					expectedDegradationRateSubset5.append(expectedDegradationRate[i])
				if len(rnaDegradationRate5_t) <= 0:
					rnaDegradationRate5.append(-1)
					expectedDegradationRateSubset5.append(-1)

				expectedDegradationRateSubset.append(expectedDegradationRate[i])

			if rnaCountsTotal[i] == 0:
				rnaDegradationRate1.append(-1)
				rnaDegradationRate2.append(-1)
				rnaDegradationRate3.append(-1)
				expectedDegradationRateSubset.append(-1)
				rnaDegradationRate4.append(-1)
				expectedDegradationRateSubset4.append(-1)
				rnaDegradationRate5.append(-1)
				expectedDegradationRateSubset5.append(-1)

		# reduction of genes
		expectedDegR = []
		predictedDegR = []
		for i in range(0,len(rnaDegradationRate3)):
			if rnaDegradationRate3[i] != -1 and rnaDegradationRate3[i] != 0 and rnaDegradationRate3[i] != 1:
				if rnaDegradationRate3[i] > 0. and rnaDegradationRate3[i] <= 0.01:
					expectedDegR.append(expectedDegradationRate[i])
					predictedDegR.append(rnaDegradationRate3[i])


		plt.figure(figsize = (8.5, 11))
		maxLine = 1.1 * max(max(expectedDegR), max(predictedDegR))

		plt.plot([0, maxLine], [0, maxLine], '--r')
		plt.plot(expectedDegR, predictedDegR, 'o', markeredgecolor = 'k', markerfacecolor = 'none')
		Correlation_ExpPred = np.corrcoef(expectedDegR, predictedDegR)[0][1]

		plt.xlabel("Expected RNA decay")
		plt.ylabel("Actual RNA decay (at final time step)")

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)


if __name__ == "__main__":
	Plot().cli()
