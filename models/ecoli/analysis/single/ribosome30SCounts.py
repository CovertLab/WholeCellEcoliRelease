"""
Plots counts of 30S rRNA, associated proteins, and complexes

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 9/5/2014
"""

from __future__ import absolute_import, division, print_function

import os

import numpy as np
from matplotlib import pyplot as plt
from six.moves import cPickle, range

from wholecell.io.tablereader import TableReader
from wholecell.utils.sparkline import sparklineAxis, setAxisMaxMinY
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from models.ecoli.analysis import singleAnalysisPlot

FONT = {
	'size':	8
	}


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Load data from KB
		sim_data = cPickle.load(open(simDataFile, "rb"))
		proteinIds = sim_data.molecule_groups.s30_proteins
		rnaIds = [sim_data.process.translation.monomer_data['rna_id'][np.where(sim_data.process.translation.monomer_data['id'] == pid)[0][0]] for pid in proteinIds]
		rRnaIds = sim_data.molecule_groups.s30_16s_rRNA
		complexIds = [sim_data.molecule_ids.s30_full_complex]

		# Load count data for mRNAs
		mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'mRNACounts'))
		mRNA_counts = mRNA_counts_reader.readColumn('mRNA_counts')
		all_mRNA_idx = {rna: i for i, rna in enumerate(mRNA_counts_reader.readAttribute('mRNA_ids'))}
		rnaIndexes = np.array([all_mRNA_idx[rna] for rna in rnaIds], np.int)
		rnaCounts = mRNA_counts[:, rnaIndexes]
		(freeProteinCounts, freeRRnaCounts, complexCounts) = read_bulk_molecule_counts(
			simOutDir, (proteinIds, rRnaIds, complexIds))
		complexCounts = complexCounts.reshape(-1, 1)

		# Load data
		main_reader = TableReader(os.path.join(simOutDir, "Main"))
		initialTime = main_reader.readAttribute("initialTime")
		time = main_reader.readColumn("time") - initialTime

		uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))

		ribosomeIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index('active_ribosome')
		activeRibosome = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, ribosomeIndex]

		plt.figure(figsize = (8.5, 15))
		plt.rc('font', **FONT)

		for idx in range(len(proteinIds)):
			rna_axis = plt.subplot(12, 3, idx + 1)

			sparklineAxis(rna_axis, time / 60., rnaCounts[:, idx], 'left', '-', 'b')
			setAxisMaxMinY(rna_axis, rnaCounts[:, idx])

			protein_axis = rna_axis.twinx()
			sparklineAxis(protein_axis, time / 60., freeProteinCounts[:, idx], 'right', '-', 'r')
			setAxisMaxMinY(protein_axis, freeProteinCounts[:, idx])

			# Component label
			rna_axis.set_xlabel(proteinIds[idx][:-3])

		for idx in range(len(rRnaIds)):
			rna_axis = plt.subplot(12, 3, idx + len(proteinIds) + 1)

			sparklineAxis(rna_axis, time / 60., freeRRnaCounts[:, idx], 'left', '-', 'b')
			setAxisMaxMinY(rna_axis, freeRRnaCounts[:, idx])

			# Component label
			rna_axis.set_xlabel(rRnaIds[idx][:-3])

		for idx in range(len(complexIds)):
			complex_axis = plt.subplot(12, 3, idx + len(proteinIds) + len(rRnaIds) + 1)

			sparklineAxis(complex_axis, time / 60., complexCounts[:, idx], 'left', '-', 'r')
			setAxisMaxMinY(complex_axis, complexCounts[:, idx])

			# Component label
			complex_axis.set_xlabel(complexIds[idx][:-3])

		# Plot number of ribosomes
		ribosome_axis = plt.subplot(12, 3, 1 + len(proteinIds) + len(rRnaIds) + len(complexIds) + 1)
		sparklineAxis(ribosome_axis, time / 60., activeRibosome, 'left', '-', 'r')
		setAxisMaxMinY(ribosome_axis, activeRibosome)
		ribosome_axis.set_xlabel('Active ribosome')

		# Save
		plt.tight_layout()

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
