"""
Plots the counts of active replisomes and their subunits over time.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 7/30/2018
"""

from __future__ import absolute_import, division, print_function

import os
from six.moves import cPickle

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from models.ecoli.analysis import singleAnalysisPlot


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)

		# Listeners used
		unique_molecule_counts_reader = TableReader(
			os.path.join(simOutDir, "UniqueMoleculeCounts")
			)
		replication_reader = TableReader(
			os.path.join(simOutDir, "ReplicationData")
			)
		main_reader = TableReader(
			os.path.join(simOutDir, "Main")
			)

		# Load counts of DNA polymerases, active replisomes, and OriC's
		unique_molecule_ids = [
			"active_replisome", 'oriC']
		unique_molecule_idx = np.array([unique_molecule_counts_reader.readAttribute(
			"uniqueMoleculeIds").index(x) for x in unique_molecule_ids])
		unique_molecule_counts = unique_molecule_counts_reader.readColumn(
			"uniqueMoleculeCounts")[:, unique_molecule_idx]

		# Load data on cell mass per origin
		criticalMassPerOriC = replication_reader.readColumn("criticalMassPerOriC")

		# Load IDs of direct replisome subunits
		replisome_subunit_ids = []
		replisome_subunit_ids.extend(
			sim_data.molecule_groups.replisome_trimer_subunits)
		replisome_subunit_ids.extend(
			sim_data.molecule_groups.replisome_monomer_subunits)

		# Load IDs of DNA polymerase III core enzyme subunits
		replisome_subunit_ids.extend(
			sim_data.process.complexation.get_monomers(
				'CPLX0-2361[c]')['subunitIds']
			)

		# Load counts of replisome subunits
		(replisome_subunit_counts,) = read_bulk_molecule_counts(simOutDir, (replisome_subunit_ids,))

		# Load time data
		initialTime = main_reader.readAttribute("initialTime")
		time = main_reader.readColumn("time") - initialTime

		# Plot figures
		fig = plt.figure()
		plot_count = 1 + len(unique_molecule_ids) + len(replisome_subunit_ids)
		fig.set_size_inches(8, 3*plot_count)
		gs = gridspec.GridSpec(plot_count, 1)

		# Plot criticalMassPerOriC
		ax = plt.subplot(gs[0, 0])
		ax.plot(time, criticalMassPerOriC)
		ax.plot(time, np.ones_like(time), color='k', linestyle='dashed')
		ax.set_xlabel("Time [s]")
		ax.set_ylabel("criticalMassPerOriC")

		# Plot counts of unique molecules
		for idx, name in enumerate(unique_molecule_ids):
			ax = plt.subplot(gs[1 + idx, 0])
			ax.plot(time, unique_molecule_counts[:, idx])
			ax.set_xlabel("Time [s]")
			ax.set_ylabel("%s\ncounts" % (name, ))

		# Plot counts of replisome subunits
		for idx, name in enumerate(replisome_subunit_ids):
			ax = plt.subplot(gs[1 + len(unique_molecule_ids) + idx, 0])
			ax.plot(time, replisome_subunit_counts[:, idx])
			ax.set_xlabel("Time [s]")
			ax.set_ylabel("%s\ncounts" % (name, ))

		fig.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close()


if __name__ == "__main__":
	Plot().cli()
