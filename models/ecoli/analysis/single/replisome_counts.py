"""
Plots the counts of active replisomes and their subunits over time.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 7/30/2018
"""

from __future__ import absolute_import
from __future__ import division

import os
import cPickle

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

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

		# Load counts of DNA polymerases, active replisomes, and OriC's
		uniqueMoleculeCounts = TableReader(
			os.path.join(simOutDir, "UniqueMoleculeCounts"))
		unique_molecule_ids = [
			"activeReplisome", "dnaPolymerase", "originOfReplication"]
		unique_molecule_idx = np.array([uniqueMoleculeCounts.readAttribute(
			"uniqueMoleculeIds").index(x) for x in unique_molecule_ids])
		unique_molecule_counts = uniqueMoleculeCounts.readColumn(
			"uniqueMoleculeCounts")[:, unique_molecule_idx]

		# Load data on cell mass per origin
		replicationData = TableReader(
			os.path.join(simOutDir, "ReplicationData"))
		criticalMassPerOriC = replicationData.readColumn("criticalMassPerOriC")

		# Load IDs of replisome subunits
		sim_data = cPickle.load(open(simDataFile, "rb"))
		replisome_subunit_ids = []
		replisome_subunit_ids.extend(
			sim_data.moleculeGroups.replisome_trimer_subunits)
		replisome_subunit_ids.extend(
			sim_data.moleculeGroups.replisome_monomer_subunits)

		# Load counts of replisome subunits
		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		moleculeIds = bulkMolecules.readAttribute("objectNames")
		replisome_subunit_idx = np.array([moleculeIds.index(x)
			for x in replisome_subunit_ids])
		replisome_subunit_counts = bulkMolecules.readColumn(
			"counts")[:, replisome_subunit_idx]

		# Load time data
		initialTime = TableReader(
			os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		time = TableReader(
			os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

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
