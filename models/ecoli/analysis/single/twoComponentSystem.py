"""
Plot two component system counts
"""

from __future__ import absolute_import, division, print_function

import os

import numpy as np
from matplotlib import pyplot as plt
from six.moves import cPickle, range

from wholecell.io.tablereader import TableReader
from wholecell.utils import units
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from models.ecoli.analysis import singleAnalysisPlot
from six.moves import zip


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):

		sim_data = cPickle.load(open(simDataFile, "rb"))
		TCS_IDS = []
		moleculeTypeOrder = ["HK", "PHOSPHO-HK", "LIGAND", "HK-LIGAND", "PHOSPHO-HK-LIGAND", "RR", "PHOSPHO-RR"]
		moleculeTypeColor = ["b", "b", "orange", "g", "g", "r", "r"]
		for system in sim_data.molecule_groups.twoComponentSystems:
			for idx, moleculeType in enumerate(moleculeTypeOrder):
				TCS_IDS.append(
					str(system["molecules"][moleculeType])
					+ sim_data.getter.get_compartment_tag(system["molecules"][moleculeType]))

		(moleculeCounts,) = read_bulk_molecule_counts(simOutDir, (TCS_IDS,))

		main_reader = TableReader(os.path.join(simOutDir, "Main"))
		initialTime = main_reader.readAttribute("initialTime")
		time = main_reader.readColumn("time") - initialTime

		# Convert molecule counts to concentrations
		nAvogadro = sim_data.constants.n_avogadro.asNumber(1 / units.mol)
		cellDensity = sim_data.constants.cell_density.asNumber(units.g / units.L)
		mass = TableReader(os.path.join(simOutDir, "Mass"))
		cellMass = (units.fg * mass.readColumn("cellMass")).asNumber(units.g)
		cellVolume = cellMass / cellDensity

		rows = 3
		cols = 5
		num_subentries = 7

		plt.figure(figsize = (8.5, 11))

		RRs = [str(x["molecules"]["RR"]) for x in sim_data.molecule_groups.twoComponentSystems]
		RRs_unique = []
		for x in RRs:
			if x not in RRs_unique:
				RRs_unique.append(x)
		RR_phosphorylation = dict(zip(RRs_unique, np.zeros(len(RRs_unique))))

		RR = []
		RRP = []
		new_RR = True


		for idx in range(len(sim_data.molecule_groups.twoComponentSystems)):
			grid_loc = idx + 1 + (cols*(num_subentries + 1))*( idx / cols)
			current_RR = str(sim_data.molecule_groups.twoComponentSystems[idx]["molecules"]["RR"])
			if RR_phosphorylation[current_RR].size == 1:
				new_RR = True

			for subentryIdx in range(len(moleculeTypeOrder)):
				if new_RR:
					if moleculeTypeOrder[subentryIdx] == "RR":
						RR[:] = moleculeCounts[:, (idx * num_subentries) + subentryIdx] / (cellVolume * nAvogadro)

					if moleculeTypeOrder[subentryIdx] == "PHOSPHO-RR":
						RRP[:] = moleculeCounts[:, (idx * num_subentries) + subentryIdx] / (cellVolume * nAvogadro)


				ax = plt.subplot(rows*(num_subentries + 2), cols, grid_loc + (cols * subentryIdx))
				ax.plot(time / 60., moleculeCounts[:, (idx * num_subentries) + subentryIdx] / (cellVolume * nAvogadro), linewidth = 1, color = moleculeTypeColor[subentryIdx])
				ax.set_title(TCS_IDS[(idx * num_subentries) + subentryIdx], fontsize = 4)

				ymin = np.amin(moleculeCounts[:, (idx * num_subentries) + subentryIdx] / (cellVolume * nAvogadro))
				ymax = np.amax(moleculeCounts[:, (idx * num_subentries) + subentryIdx] / (cellVolume * nAvogadro))
				self.set_ylim(ax, ymin, ymax)
				ax.set_yticks([ymin, ymax])
				ax.set_yticklabels(["%0.2e" % ymin, "%0.2e" % ymax])
				ax.spines['top'].set_visible(False)
				ax.spines['bottom'].set_visible(False)
				ax.xaxis.set_ticks_position('none')
				ax.tick_params(which = 'both', direction = 'out', labelsize = 4)
				ax.set_xlim(left = -2)
				ax.set_xticks([])

			# Monitor percentage of TF phosphorylation
			if new_RR:

				RR_phosphorylation[current_RR] = np.array(RRP) / (np.array(RRP) + np.array(RR))
				new_RR = False


		# Plot percentage of TF phosphorylation
		for idx_TF, idx_subentry in enumerate(np.arange(100, 135, 5)): # Using the last sub-entries
			phosphorylation = RR_phosphorylation[RRs_unique[idx_TF]]
			ax = plt.subplot(rows*(num_subentries + 2), cols, idx_subentry)
			ax.plot(time / 60., phosphorylation, linewidth = 1, color = "grey")
			ax.set_title(RRs_unique[idx_TF], fontsize = 4)

			if np.any(np.isfinite(phosphorylation)):
				ymin = np.nanmin(phosphorylation)
				ymax = np.nanmax(phosphorylation)
				if ymin == ymax:
					ymin, ymax = -0.001, 0.001
			else:
				ymin = 0
				ymax = 1
			self.set_ylim(ax, ymin, ymax)
			ax.set_yticks([ymin, ymax])
			ax.set_yticklabels(["%0.2e" % ymin, "%0.2e" % ymax])
			ax.spines['top'].set_visible(False)
			ax.spines['bottom'].set_visible(False)
			ax.xaxis.set_ticks_position('none')
			ax.tick_params(which = 'both', direction = 'out', labelsize = 4)
			ax.set_xlim(left = -2)
			ax.set_xticks([])

		# phosphorylation percentage section
		ax = plt.subplot(rows*(num_subentries + 2), cols, 95)
		ax.set_title("Fraction of response regulator phosphorylation", fontsize = 6)
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.spines['left'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.xaxis.set_ticks_position('none')
		ax.yaxis.set_ticks_position('none')
		ax.set_yticks([])
		ax.set_yticklabels([])
		ax.set_xticks([])
		ax.set_xticklabels([])

		plt.subplots_adjust(hspace = 1, wspace = 1)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
