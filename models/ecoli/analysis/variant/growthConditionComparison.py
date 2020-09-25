"""
Plot to compare cell properties across different growth conditions similar to Dennis and Bremer. 1996. Fig 2
Multiple generations of the same variant will be plotted together
"""

from __future__ import absolute_import, division, print_function


import os

import numpy as np
from matplotlib import pyplot as plt
from six.moves import cPickle

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import variantAnalysisPlot

NT_MW = 487.0
PROTEIN_MW = 110.0


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		ap = AnalysisPaths(inputDir, variant_plot = True)
		all_cells = ap.get_cells()

		rnaToProteinDict = {}
		dnaToProteinDict = {}
		elngRateDict = {}
		stableRnaFractionDict = {}
		doublingPerHourDict = {}

		variantSimDataFile = ap.get_variant_kb(ap.get_variants()[0])
		sim_data = cPickle.load(open(variantSimDataFile, "rb"))
		nAvogadro = sim_data.constants.n_avogadro.asNumber()
		chromMass = (sim_data.getter.get_mass(sim_data.molecule_ids.full_chromosome) / sim_data.constants.n_avogadro).asNumber()

		for simDir in all_cells:
			simOutDir = os.path.join(simDir, "simOut")
			variant = int(simDir[simDir.rfind('generation_')-14:simDir.rfind('generation_')-8])

			mass = TableReader(os.path.join(simOutDir, "Mass"))

			protein = mass.readColumn("proteinMass") * 10**-15
			rna = mass.readColumn("rnaMass") * 10**-15
			dna = mass.readColumn("dnaMass") * 10**-15

			growthRate = mass.readColumn("instantaniousGrowthRate")
			doublingTime = np.nanmean(np.log(2) / growthRate / 60)

			rnaNT = rna / NT_MW * nAvogadro
			proteinAA = protein / PROTEIN_MW * nAvogadro

			# Count chromosome equivalents
			chromEquivalents = dna / chromMass

			# Load ribosome data
			ribosomeDataFile = TableReader(os.path.join(simOutDir, "RibosomeData"))
			actualElongations = ribosomeDataFile.readColumn("actualElongations")
			ribosomeDataFile.close()

			transcriptDataFile = TableReader(os.path.join(simOutDir, "TranscriptElongationListener"))
			rnaSynth = transcriptDataFile.readColumn("countRnaSynthesized")
			isTRna = sim_data.process.transcription.rna_data['is_tRNA']
			isRRna = sim_data.process.transcription.rna_data['is_rRNA']
			stableRnaSynth = np.sum(rnaSynth[:,isTRna], axis=1) + np.sum(rnaSynth[:,isRRna], axis=1)
			totalRnaSynth = np.sum(rnaSynth, axis=1).astype(float)
			rnaFraction = stableRnaSynth / totalRnaSynth

			uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))

			ribosomeIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index('active_ribosome')
			activeRibosome = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, ribosomeIndex]

			uniqueMoleculeCounts.close()

			timeStepSec = TableReader(os.path.join(simOutDir, "Main")).readColumn("timeStepSec")

			if variant in rnaToProteinDict:
				rnaToProteinDict[variant] = np.append(rnaToProteinDict[variant], rnaNT / (proteinAA / 100))
				dnaToProteinDict[variant] = np.append(dnaToProteinDict[variant], chromEquivalents / (proteinAA / 10**9))
				elngRateDict[variant] = np.append(elngRateDict[variant], (actualElongations / activeRibosome / timeStepSec)[3:])
				stableRnaFractionDict[variant] = np.append(stableRnaFractionDict[variant], np.asarray(rnaFraction)[~np.isnan(rnaFraction)])
				doublingPerHourDict[variant] = np.append(doublingPerHourDict[variant], 60 / doublingTime)
			else:
				rnaToProteinDict[variant] = rnaNT / (proteinAA / 100)
				dnaToProteinDict[variant] = chromEquivalents / (proteinAA / 10**9)
				elngRateDict[variant] = (actualElongations / activeRibosome / timeStepSec)[3:]
				stableRnaFractionDict[variant] = np.asarray(rnaFraction)[~np.isnan(rnaFraction)]
				doublingPerHourDict[variant] = 60 / doublingTime

		rnaToProtein = []
		dnaToProtein = []
		elngRate = []
		stableRnaFraction = []
		doublingPerHour = []

		for key in rnaToProteinDict:
			rnaToProtein += [rnaToProteinDict[key]]
			dnaToProtein += [dnaToProteinDict[key]]
			elngRate += [elngRateDict[key]]
			stableRnaFraction += [stableRnaFractionDict[key]]
			doublingPerHour += [np.mean(doublingPerHourDict[key])]

		plt.figure(figsize = (8.5, 11))

		sp = plt.subplot(4,1,1)
		sp.violinplot(rnaToProtein, positions=doublingPerHour, showmeans=True)
		sp.set_ylabel("RNA to Protein\n(nuc/100 aa)")

		sp = plt.subplot(4,1,2)
		sp.violinplot(dnaToProtein, positions=doublingPerHour, showmeans=True)
		sp.set_ylabel("DNA to Protein\n(chrom eq/10^9 aa)")

		sp = plt.subplot(4,1,3)
		sp.violinplot(elngRate, positions=doublingPerHour, showmeans=True)
		sp.set_ylabel("Ribosome Elongation\nRate (aa/s)")

		sp = plt.subplot(4,1,4)
		sp.violinplot(stableRnaFraction, positions=doublingPerHour, showmeans=True)
		sp.set_ylabel("Rate Stable RNA to\nRate Total RNA")
		sp.set_xlabel("Doublings per Hour")

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
