"""
Plot rna synthesis probabilities

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 9/9/2016
"""

from __future__ import absolute_import
from __future__ import division

import cPickle
import os
from itertools import izip

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import filepath


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, 'simOutDir does not currently exist as a directory'

		filepath.makedirs(plotOutDir)

		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)

		# Load info from sim_data
		isMRna = sim_data.process.transcription.rnaData["isMRna"]
		isRRna = sim_data.process.transcription.rnaData["isRRna"]
		isTRna = sim_data.process.transcription.rnaData["isTRna"]
		nutrients = sim_data.conditions[sim_data.condition]["nutrients"]
		synth_prob_fractions = sim_data.process.transcription.rnaSynthProbFraction[nutrients]

		# Get "average" synthesis probability fractions set by the parca
		mrna_avg_synth_prob = synth_prob_fractions["mRna"]
		trna_avg_synth_prob = synth_prob_fractions["tRna"]
		rrna_avg_synth_prob = synth_prob_fractions["rRna"]

		# Listeners used
		rna_synth_prob_reader = TableReader(os.path.join(simOutDir, "RnaSynthProb"))
		mass_reader = TableReader(os.path.join(simOutDir, "Mass"))

		# Load data
		rna_synth_prob = rna_synth_prob_reader.readColumn('rnaSynthProb')
		time = rna_synth_prob_reader.readColumn('time')
		mrna_synth_prob = rna_synth_prob[:, isMRna].sum(axis = 1)
		trna_synth_prob = rna_synth_prob[:, isTRna].sum(axis = 1)
		rrna_synth_prob = rna_synth_prob[:, isRRna].sum(axis = 1)

		n_mrnas = np.sum(isMRna)
		n_trnas = np.sum(isTRna)
		n_rrnas = np.sum(isRRna)

		mrna_mass = mass_reader.readColumn("mRnaMass")
		trna_mass = mass_reader.readColumn("tRnaMass")
		rrna_mass = mass_reader.readColumn("rRnaMass")

		# Plot
		fig = plt.figure(figsize = (10, 15))

		synth_probs = [mrna_synth_prob, trna_synth_prob, rrna_synth_prob]
		avg_synth_probs = [
			mrna_avg_synth_prob,
			trna_avg_synth_prob,
			rrna_avg_synth_prob,
			]
		normalized_mass = [
			mrna_mass/mrna_mass[0],
			trna_mass/trna_mass[0],
			rrna_mass/rrna_mass[0],
			]
		subplot_titles = [
			"mRNA\n(sum of %s mRNAs)" % n_mrnas,
			"tRNA\n(sum of %s tRNAs)" % n_trnas,
			"rRNA\n(sum of %s rRNAs)" % n_rrnas,
			]

		for index, (title, synth_prob, avg_synth_prob, mass) in enumerate(izip(
				subplot_titles, synth_probs, avg_synth_probs, normalized_mass)):

			ax1 = plt.subplot(3, 1, index + 1)
			ax1.plot(time[1:], synth_prob[1:], label="Sum of synthesis probabilities")
			ax1.axhline(avg_synth_prob,
				linestyle="--", color='k', linewidth=3,
				label="Fit fraction of transcription probabilities")

			ax1.set_title(title)
			ax1.set_xlim([time[0], time[-1]])
			ax1.set_ylim([0, synth_prob.max()*1.2])
			ax1.set_xlabel("Time [s]")
			ax1.set_ylabel("Transcription probability")
			ax1.legend(loc=2)

			ax2 = ax1.twinx()
			ax2.plot(time[1:], mass[1:], label="Normalized mass", color='r')
			ax2.set_ylabel("Mass (normalized by t=0)")
			ax2.legend(loc=1)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
