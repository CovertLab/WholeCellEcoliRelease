"""
Plot trp regulation

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/17/2016
"""

from __future__ import absolute_import, division, print_function

import os

import numpy as np
from matplotlib import pyplot as plt
import six
from six.moves import cPickle, range

from wholecell.io.tablereader import TableReader
from wholecell.utils import units
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from models.ecoli.analysis import singleAnalysisPlot


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Load data from KB
		sim_data = cPickle.load(open(simDataFile, "rb"))
		nAvogadro = sim_data.constants.n_Avogadro
		cellDensity = sim_data.constants.cellDensity

		# Load time
		main_reader = TableReader(os.path.join(simOutDir, "Main"))
		initialTime = main_reader.readAttribute("initialTime")
		time = main_reader.readColumn("time") - initialTime

		# Load mass data
		# Total cell mass is needed to compute concentrations (since we have cell density)
		# Protein mass is needed to compute the mass fraction of the proteome that is trpA
		massReader = TableReader(os.path.join(simOutDir, "Mass"))

		cellMass = units.fg * massReader.readColumn("cellMass")
		proteinMass = units.fg * massReader.readColumn("proteinMass")

		# Load data from RnaSynthProb listener
		rna_synth_prob_reader = TableReader(
			os.path.join(simOutDir, "RnaSynthProb"))
		tf_ids = rna_synth_prob_reader.readAttribute("tf_ids")
		rna_idx = {rna: i for i, rna in enumerate(rna_synth_prob_reader.readAttribute("rnaIds"))}
		n_bound_TF_per_TU = rna_synth_prob_reader.readColumn(
			"n_bound_TF_per_TU").reshape(
			(-1, len(rna_idx), len(tf_ids)))

		# Get indexes of trpR and its target RNAs
		trpRIndex = tf_ids.index("CPLX-125")
		target_ids = six.viewkeys(sim_data.tf_to_fold_change["CPLX-125"])
		target_idx = np.array(
			[rna_idx[target_id + "[c]"] for target_id in target_ids])

		trpAProteinId = ["TRYPSYN-APROTEIN[c]"]
		bulk_ids = (
			["TRP[c]"],  # intracellular trp
			["CPLX-125[c]"],  # active trpR (that isn't promoter bound)
			["PC00007[c]"],  # inactive trpR
			["PD00423[c]"],  # monomeric trpR
			trpAProteinId,  # monomeric trpA
			["TRYPSYN[c]"],  # complexed trpA
			)
		(trpCounts, trpRActiveCounts, trpRInactiveCounts, trpRMonomerCounts,
			trpAProteinCounts, trpABComplexCounts
			) = read_bulk_molecule_counts(simOutDir, bulk_ids)

		trpMols = 1. / nAvogadro * trpCounts
		volume = cellMass / cellDensity
		trpConcentration = trpMols * 1. / volume

		# Get the promoter-bound status for all regulated genes
		tfBoundCounts = n_bound_TF_per_TU[:, target_idx, trpRIndex]

		# Compute total counts of trpA in monomeric and complexed form
		# (we know the stoichiometry)
		trpAProteinTotalCounts = trpAProteinCounts + 2 * trpABComplexCounts

		# Compute the trpA mass in the cell
		trpAMw = sim_data.getter.get_mass(trpAProteinId)
		trpAMass = 1. / nAvogadro * trpAProteinTotalCounts * trpAMw

		# Compute the proteome mass fraction
		proteomeMassFraction = trpAMass.asNumber(units.fg) / proteinMass.asNumber(units.fg)

		# Get the synthesis probability for all regulated genes
		synthProbIds = [target + "[c]" for target in sim_data.tf_to_fold_change["CPLX-125"].keys()]
		synthProbIndex = np.array([rna_idx[x] for x in synthProbIds])
		synthProbs = rna_synth_prob_reader.readColumn("rnaSynthProb")[:, synthProbIndex]

		tf_ids = rna_synth_prob_reader.readAttribute("tf_ids")
		trpRIndex = [i for i, tf in enumerate(tf_ids) if tf == "CPLX-125"][0]
		trpRBound = rna_synth_prob_reader.readColumn("nActualBound")[:,trpRIndex]

		# Calculate total trpR - active, inactive, bound and monomeric
		trpRTotalCounts = 2 * (trpRActiveCounts + trpRInactiveCounts + trpRBound) + trpRMonomerCounts

		# Compute moving averages
		width = 100

		tfBoundCountsMA = np.array([np.convolve(tfBoundCounts[:,i], np.ones(width) / width, mode = "same")
				for i in range(tfBoundCounts.shape[1])]).T
		synthProbsMA = np.array([np.convolve(synthProbs[:,i], np.ones(width) / width, mode = "same")
				for i in range(synthProbs.shape[1])]).T

		plt.figure(figsize = (8.5, 11))

		##############################################################
		ax = plt.subplot(6, 1, 1)
		ax.plot(time, trpConcentration.asNumber(units.umol / units.L))
		plt.ylabel("Internal TRP Conc. [uM]", fontsize = 6)

		ymin = np.amin(trpConcentration.asNumber(units.umol / units.L) * 0.9)
		ymax = np.amax(trpConcentration.asNumber(units.umol / units.L) * 1.1)
		if ymin != ymax:
			ax.set_ylim([ymin, ymax])
			ax.set_yticks([ymin, ymax])
			ax.set_yticklabels(["%0.0f" % ymin, "%0.0f" % ymax])
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.xaxis.set_ticks_position('none')
		ax.tick_params(which = 'both', direction = 'out', labelsize = 6)
		ax.set_xticks([])
		##############################################################

		##############################################################
		ax = plt.subplot(6, 1, 2)
		ax.plot(time, trpRActiveCounts)
		ax.plot(time, trpRInactiveCounts)
		ax.plot(time, trpRTotalCounts)
		plt.ylabel("TrpR Counts", fontsize = 6)
		plt.legend(["Active (dimer)", "Inactive (dimer)", "Total (monomeric)"], fontsize = 6)

		ymin = min(np.amin(trpRActiveCounts * 0.9), np.amin(trpRInactiveCounts * 0.9))
		ymax = np.amax(trpRTotalCounts * 1.1)
		if ymin != ymax:
			ax.set_ylim([ymin, ymax])
			ax.set_yticks([ymin, ymax])
			ax.set_yticklabels(["%0.0f" % ymin, "%0.0f" % ymax])
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.xaxis.set_ticks_position('none')
		ax.tick_params(which = 'both', direction = 'out', labelsize=6)
		ax.set_xticks([])
		##############################################################

		##############################################################
		ax = plt.subplot(6, 1, 3)
		ax.plot(time, tfBoundCountsMA)
		plt.ylabel("TrpR Bound To Promoters\n(Moving Average)", fontsize = 6)

		ymin = np.amin(tfBoundCountsMA * 1.)
		ymax = np.amax(tfBoundCountsMA * 1.)
		if ymin != ymax:
			ax.set_ylim([ymin, ymax])
			ax.set_yticks([ymin, ymax])
			ax.set_yticklabels(["%0.0f" % ymin, "%0.0f" % ymax])
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.xaxis.set_ticks_position('none')
		ax.tick_params(which = 'both', direction = 'out', labelsize = 6)
		ax.set_xticks([])
		##############################################################

		##############################################################
		ax = plt.subplot(6, 1, 4)
		ax.plot(time, synthProbsMA)
		plt.ylabel("Regulated Gene Synthesis Prob.\n(Moving Average)", fontsize = 6)

		ymin = np.amin(synthProbsMA[1:] * 0.9)
		ymax = np.amax(synthProbsMA[1:] * 1.1)
		if ymin != ymax:
			ax.set_ylim([ymin, ymax])
			ax.set_yticks([ymin, ymax])
			ax.set_yticklabels(["%0.2e" % ymin, "%0.2e" % ymax])
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.xaxis.set_ticks_position('none')
		ax.tick_params(which = 'both', direction = 'out', labelsize = 6)
		ax.set_xticks([])
		##############################################################

		##############################################################
		ax = plt.subplot(6, 1, 5)
		ax.plot(time, trpAProteinTotalCounts)
		plt.ylabel("TrpA Counts", fontsize = 6)

		ymin = np.amin(trpAProteinTotalCounts * 0.9)
		ymax = np.amax(trpAProteinTotalCounts * 1.1)
		if ymin != ymax:
			ax.set_ylim([ymin, ymax])
			ax.set_yticks([ymin, ymax])
			ax.set_yticklabels(["%0.0f" % ymin, "%0.0f" % ymax])
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.xaxis.set_ticks_position('none')
		ax.tick_params(which = 'both', direction = 'out', labelsize = 6)
		ax.set_xticks([])
		##############################################################


		##############################################################
		ax = plt.subplot(6, 1, 6)
		ax.plot(time, proteomeMassFraction)
		plt.ylabel("TrpA Mass Fraction of Proteome", fontsize = 6)

		ymin = np.amin(proteomeMassFraction * 0.9)
		ymax = np.amax(proteomeMassFraction * 1.1)
		if ymin != ymax:
			ax.set_ylim([ymin, ymax])
			ax.set_yticks([ymin, ymax])
			ax.set_yticklabels(["%0.2e" % ymin, "%0.2e" % ymax])
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.xaxis.set_ticks_position('none')
		ax.tick_params(which = 'both', direction = 'out', labelsize = 6)
		ax.set_xticks([])
		##############################################################

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
