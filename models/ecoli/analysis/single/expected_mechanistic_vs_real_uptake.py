import os
from six.moves import cPickle

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):

		# Amino acid IDs
		sim_data = cPickle.load(open(simDataFile, "rb"))
		aa_ids = sim_data.molecule_groups.amino_acids

		# Amino acid exchanges fluxes
		main_reader = TableReader(os.path.join(simOutDir, "Main"))
		initial_time = main_reader.readAttribute("initialTime")
		time = main_reader.readColumn("time") - initial_time

		# Uptake fluxes
		fba_results = TableReader(os.path.join(simOutDir, "FBAResults"))
		external_exchange_fluxes = fba_results.readColumn("externalExchangeFluxes")
		external_molecule_ids = fba_results.readAttribute("externalMoleculeIDs")

		#Transporter counts
		monomers = TableReader(os.path.join(simOutDir,"BulkMolecules"))
		monomer_names = monomers.readAttribute('objectNames')
		monomer_counts = monomers.readColumn('counts')

		enzyme_kinetics_reader = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))
		counts_to_molar = enzyme_kinetics_reader.readColumn('countsToMolar')

		# kcats to aa mapping
		kcats = {aa: kcat for aa, kcat in zip(aa_ids, sim_data.process.metabolism.import_kcats_per_aa)}

		# Cell Dry mass at each time step
		mass = TableReader(os.path.join(simOutDir, "Mass"))
		dry_mass = mass.readColumn('dryMass')
		cell_volume = mass.readColumn('cellVolume')

		# Map aas to their transporters
		aa_to_importers = sim_data.process.metabolism.aa_to_importers

		# Plot 1 - Expected vs real uptakes
		rows = 6
		cols = 4
		fig = plt.figure(figsize=(8, 11.5))

		for plot_index, aa in enumerate(aa_ids):
			ax = plt.subplot(rows, cols, plot_index + 1)

			# Get actual fluxes
			old_aa_tag = aa
			if not aa.startswith("L-SELENOCYSTEINE"):
				aa = aa[: -3] + "[p]"
			if aa in external_molecule_ids:
				aa_flux = external_exchange_fluxes[:, external_molecule_ids.index(aa)]
			else:
				aa_flux = np.zeros(len(time))

			#calculate expected rate based on kcats and transporters
			rate_= kcats[old_aa_tag] * np.ones(len(aa_flux))
			t_counts = np.zeros(len(aa_flux))
			for i in set(aa_to_importers[old_aa_tag]):
				t_counts += monomer_counts[:, monomer_names.index(i)]

			# Right now, code breaks because CYS has 0 transporters
			rate_ *= t_counts

			#rate in mol/g/hr
			rate_*= (-3600 * 1000) / (dry_mass * 1e-15 * sim_data.constants.n_avogadro.asNumber())

			#Plot, orange is target flux and blue is actual flux
			ax.plot(time / 60., aa_flux, linewidth = 1, label = 'Actual Flux')
			ax.plot(time / 60., rate_,  linewidth = 1, label = 'Expected Flux')
			ax.set_xlabel("Time (min)", fontsize = 6)
			ax.set_ylabel("mmol/gDCW/hr", fontsize = 6)
			ax.set_title("%s" % aa, fontsize = 6, y = 1.1)
			ax.tick_params(which = "both", direction = "out", labelsize = 6)

		plt.rc("font", size = 6)
		plt.suptitle("External exchange fluxes of amino acids", fontsize = 10)

		plt.subplots_adjust(hspace = 1, wspace = 1)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")

		# Plot 2 - Transporters counts
		rows = 8
		cols = 6
		fig = plt.figure(figsize = (8, 11.5))

		mmIDs = sim_data.process.metabolism.aa_importer_names

		for plot_index, m in enumerate(mmIDs):
			ax = plt.subplot(rows, cols, plot_index + 1)

			if m in monomer_names:
				mC = monomer_counts[:, monomer_names.index(m)]
			else:
				mC = np.zeros(len(time))

			ax.plot(time / 60., mC * counts_to_molar)
			ax.set_xlabel("Time (min)", fontsize = 6)
			ax.set_ylabel("Molar", fontsize = 6)
			ax.set_title("%s" % m, fontsize = 6, y = 1.1)
			ax.tick_params(which = "both", direction = "out", labelsize = 6)

		plt.rc("font", size = 6)
		plt.suptitle("Counts of monomers of interest", fontsize = 10)

		plt.subplots_adjust(hspace = 1, wspace = 1)
		exportFigure(plt, plotOutDir, 'transporter_conc', metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
