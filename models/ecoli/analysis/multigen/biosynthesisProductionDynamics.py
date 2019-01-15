"""
Plot the protein and RNA levels, as well as synthesis probabilities,
for regulated genes that are involved in AA biosynthesis

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/17/2016
"""

from __future__ import absolute_import
from __future__ import division

import os

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(seedOutDir):
			raise Exception, "seedOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)

		allDirs = ap.get_cells()

		tfs = [
			"putA", "aldA", "gdhA", "carA", "carB", "argD",
			"argE", "aroG", "aroF", "tyrA", "tyrB", "trpA",
			"trpB", "trpC", "trpD", "trpE", "aroH"
			]

		tfToRNAId = {
			"putA":	"EG10801_RNA[c]",
			"aldA":	"EG10035_RNA[c]",
			"gdhA":	"EG10372_RNA[c]",
			"carA":	"EG10134_RNA[c]",
			"carB":	"EG10135_RNA[c]",
			"argD":	"EG10066_RNA[c]",
			"argE":	"EG11286_RNA[c]",
			"aroG":	"EG10079_RNA[c]",
			"aroF":	"EG10078_RNA[c]",
			"tyrA":	"EG11039_RNA[c]",
			"tyrB":	"EG11040_RNA[c]",
			"trpA":	"EG11024_RNA[c]",
			"trpB":	"EG11025_RNA[c]",
			"trpC":	"EG11026_RNA[c]",
			"trpD":	"EG11027_RNA[c]",
			"trpE":	"EG11028_RNA[c]",
			"aroH":	"EG10080_RNA[c]",
		}

		tfToMonomerId = {
			"putA":	"PUTA-MONOMER[c]",
			"aldA":	"LACTALDDEHYDROG-MONOMER[c]",
			"gdhA":	"GDHA-MONOMER[c]",
			"carA":	"CARBPSYN-SMALL[c]",
			"carB":	"CARBPSYN-LARGE[c]",
			"argD":	"ACETYLORNTRANSAM-MONOMER[c]",
			"argE":	"ACETYLORNDEACET-MONOMER[c]",
			"aroG":	"AROG-MONOMER[c]",
			"aroF":	"AROF-MONOMER[c]",
			"tyrA":	"CHORISMUTPREPHENDEHYDROG-MONOMER[c]",
			"tyrB":	"TYRB-MONOMER[c]",
			"trpA":	"TRYPSYN-APROTEIN[c]",
			"trpB":	"TRYPSYN-BPROTEIN[c]",
			"trpC":	"PRAI-IGPS[c]",
			"trpD":	"ANTHRANSYNCOMPII-MONOMER[c]",
			"trpE":	"ANTHRANSYNCOMPI-MONOMER[c]",
			"aroH":	"AROH-MONOMER[c]",
		}

		tfToComplexId = {
			"putA":	"PUTA-CPLX[c]",
			"aldA":	"ALD-CPLX[c]",
			"gdhA":	"GDHA-CPLX[c]",
			"carA":	"CARBPSYN-CPLX[c]",
			"carB":	"CARBPSYN-CPLX[c]",
			"argD":	"ACETYLORNTRANSAM-CPLX[c]",
			"argE":	"ACETYLORNDEACET-CPLX[c]",
			"aroG":	"AROG-CPLX[c]",
			"aroF":	"AROF-CPLX[c]",
			"tyrA":	"CHORISMUTPREPHENDEHYDROG-CPLX[c]",
			"tyrB":	"TYRB-DIMER[c]",
			"trpA":	"TRYPSYN[c]",
			"trpB":	"TRYPSYN[c]",
			"trpD":	"ANTHRANSYN-CPLX[c]",
			"trpE":	"ANTHRANSYN-CPLX[c]",
			"aroH":	"AROH-CPLX[c]",
		}

		tfToComplexStoich = {
			"putA":	2,
			"aldA":	4,
			"gdhA":	6,
			"carA":	2,
			"carB":	2,
			"argD":	2,
			"argE":	2,
			"aroG":	4,
			"aroF":	2,
			"tyrA":	2,
			"tyrB":	2,
			"trpA":	2,
			"trpB":	2,
			"trpD":	2,
			"trpE":	2,
			"aroH":	2,
		}

		nTfs = len(tfs)


		plt.figure(figsize = (8.5, 11))

		for simDir in allDirs:
			simOutDir = os.path.join(simDir, "simOut")
			# Load time
			time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")

			# Load data from bulk molecules
			bulkMoleculesReader = TableReader(os.path.join(simOutDir, "BulkMolecules"))
			bulkMoleculeIds = bulkMoleculesReader.readAttribute("objectNames")
			bulkMoleculeCounts = bulkMoleculesReader.readColumn("counts")
			bulkMoleculesReader.close()

			# Get the synthesis probability for all regulated genes
			rnaSynthProbReader = TableReader(os.path.join(simOutDir, "RnaSynthProb"))
			rnaSynthProbIds = rnaSynthProbReader.readAttribute("rnaIds")
			synthProbs = rnaSynthProbReader.readColumn("rnaSynthProb")
			rnaSynthProbReader.close()

			for tfIdx, tf in enumerate(tfs):
				monomerId = tfToMonomerId[tf]
				monomerIdx = bulkMoleculeIds.index(monomerId)
				monomerCounts = bulkMoleculeCounts[:, monomerIdx].copy()
				if tf in tfToComplexId:
					complexId = tfToComplexId[tf]
					complexIdx = bulkMoleculeIds.index(complexId)
					complexCounts = bulkMoleculeCounts[:, complexIdx].copy()

					monomerCounts += (tfToComplexStoich[tf] * complexCounts)

				rnaId = tfToRNAId[tf]
				rnaIdx = bulkMoleculeIds.index(rnaId)
				rnaCounts = bulkMoleculeCounts[:, rnaIdx].copy()

				synthProbIdx = rnaSynthProbIds.index(rnaId)
				synthProb = synthProbs[:, synthProbIdx].copy()

				# Compute moving averages
				width = 100

				synthProbMA = np.convolve(synthProb, np.ones(width) / width, mode = "same")

				##############################################################
				ax = self.subplot(nTfs, 3, tfIdx * 3 + 1)
				ax.plot(time, monomerCounts, color = "b")
				plt.title("%s counts" % monomerId, fontsize = 8)

				ymin, ymax = ax.get_ylim()
				ax.set_yticks([ymin, ymax])
				ax.set_yticklabels(["%0.0f" % ymin, "%0.0f" % ymax])
				ax.spines['top'].set_visible(False)
				ax.spines['bottom'].set_visible(False)
				ax.xaxis.set_ticks_position('none')
				ax.tick_params(which = 'both', direction = 'out', labelsize = 8)
				ax.set_xticks([])
				##############################################################

				##############################################################
				ax = self.subplot(nTfs, 3, tfIdx * 3 + 2)
				ax.plot(time, rnaCounts, color = "b")
				plt.title("%s counts" % rnaId, fontsize = 8)

				ymin, ymax = ax.get_ylim()
				ax.set_yticks([ymin, ymax])
				ax.set_yticklabels(["%0.0f" % ymin, "%0.0f" % ymax])
				ax.spines['top'].set_visible(False)
				ax.spines['bottom'].set_visible(False)
				ax.xaxis.set_ticks_position('none')
				ax.tick_params(which = 'both', direction = 'out', labelsize = 8)
				ax.set_xticks([])
				##############################################################

				##############################################################
				ax = self.subplot(nTfs, 3, tfIdx * 3 + 3)
				ax.plot(time, synthProb, color = "b")
				ax.plot(time, synthProbMA, color = "k")
				plt.title("%s synth prob" % rnaId, fontsize = 8)

				ymin, ymax = ax.get_ylim()
				ax.set_yticks([ymin, ymax])
				ax.set_yticklabels(["%0.2e" % ymin, "%0.2e" % ymax])
				ax.spines['top'].set_visible(False)
				ax.spines['bottom'].set_visible(False)
				ax.xaxis.set_ticks_position('none')
				ax.tick_params(which = 'both', direction = 'out', labelsize = 8)
				ax.set_xticks([])
				##############################################################

		plt.tight_layout()

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
