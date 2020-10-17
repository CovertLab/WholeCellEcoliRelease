from __future__ import absolute_import, division, print_function

from six.moves import cPickle
import os

from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot
from six.moves import zip


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		return


		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = cPickle.load(f)

		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)

		# Get all cells
		allDir = ap.get_cells()

		# Constants
		aa_c_ids = sim_data.molecule_groups.amino_acids
		aa_p_ids = [aa_id.replace("[c]", "[p]") for aa_id in aa_c_ids]


		for simDir in allDir:
			simOutDir = os.path.join(simDir, "simOut")




			# import ipdb; ipdb.set_trace()

			# Listeners used
			main_reader = TableReader(os.path.join(simOutDir, 'Main'))
			transport_reader = TableReader(os.path.join(simOutDir, 'TransportListener'))
			fba_reader = TableReader(os.path.join(simOutDir, "FBAResults"))
			mass_reader = TableReader(os.path.join(simOutDir, "Mass"))

			# Load time
			initial_time = main_reader.readAttribute('initialTime')
			time = main_reader.readColumn('time')[START_TIME_STEP:] - initial_time
			time_step = main_reader.readColumn('timeStepSec')[START_TIME_STEP:]

			# growth rate
			expected_doubling_time = sim_data.doubling_time
			massNames = ['instantaniousGrowthRate']
			growthRate = mass_reader.readColumn(massNames[0])[START_TIME_STEP:]
			growthRate = (1 / units.s) * growthRate
			doublingTime = 1 / growthRate * np.log(2)
			avgDoublingTime = doublingTime[1:].asNumber(units.min).mean()
			stdDoublingTime = doublingTime[1:].asNumber(units.min).std()
			doublingTime = np.array(doublingTime.asNumber())

			# transport listener
			transport_reaction_ids = transport_reader.readAttribute('reaction_ids')
			reactions = transport_reader.readAttribute('reactions')

			molecule_ids = transport_reader.readAttribute('molecule_ids')
			internal_molecule_ids = transport_reader.readAttribute('internal_molecule_ids')
			delta_molecules = transport_reader.readColumn('delta_molecules')[START_TIME_STEP:, :]
			exchange_fluxes = transport_reader.readColumn('exchange_fluxes')[START_TIME_STEP:, :]
			reaction_fluxes = transport_reader.readColumn('reaction_fluxes')[START_TIME_STEP:, :]
			molecule_concentrations = transport_reader.readColumn('molecule_concentrations')[START_TIME_STEP:, :]

			delta_molecules_dict = dict(zip(internal_molecule_ids, delta_molecules.T))
			exchange_fluxes_dict = dict(zip(molecule_ids, exchange_fluxes.T))
			concentrations_dict = dict(zip(molecule_ids, molecule_concentrations.T))
			kinetic_flux_dict = dict(zip(transport_reaction_ids, reaction_fluxes.T))

			# fba listener
			exFlux = fba_reader.readColumn("externalExchangeFluxes")
			exMolec = fba_reader.readAttribute("externalMoleculeIDs")
			glc_flux = -1. * exFlux[START_TIME_STEP:, exMolec.index("GLC[p]")]

			reaction_ids = fba_reader.readAttribute("reactionIDs")
			reaction_fluxes = fba_reader.readColumn("reactionFluxes")[START_TIME_STEP:,
							  :] * 1e-3  # convert to mol/L/s from mmol/L/s # units.mmol / units.L / units.s
			fba_flux_dict = dict(zip(reaction_ids, reaction_fluxes.T))

			# make dictionary with transporters, substrates ids for all reactions.
			transport_reaction_dict = {reaction_id: {} for reaction_id in transport_reaction_ids}
			transporters_list = []
			for reaction in reactions:
				reaction_id = reaction['reaction id']
				stoich = reaction['stoichiometry']

				transporters = []
				if 'catalyzed by' in reaction:
					transporters_no_compartment = reaction['catalyzed_by']

					# get transporter's compartment
					transporters = []
					for mol_no_compart in transporters_no_compartment:
						transporter = [mol_id for mol_id in molecule_ids if mol_no_compart in mol_id]
						if transporter:
							transporters.append(transporter[0])
							transporters_list.append(transporter[0])

				transport_reaction_dict[reaction_id] = {
					'stoichiometry': stoich,
					'transporters': transporters,
				}

			transporters_list = list(set(transporters_list))

			# get transporter concentrations from concentrations_dict
			transporters_concs = np.zeros((len(transporters_list), len(time)))
			for index, transporter in enumerate(transporters_list):
				transporters_concs[index, :] = concentrations_dict[transporter]

			transporters_concs_avg = transporters_concs.mean(axis=0)

			# TODO -- get AA exchange flux -- either from fba or from transport...
			aa_fluxes = np.zeros((len(aa_p_ids), len(time)))
			for index, aa_id in enumerate(aa_p_ids):
				aa_fluxes[index, :] = exchange_fluxes_dict[aa_id]

			aa_fluxes_avg = aa_fluxes.mean(axis=0)
			aa_fluxes_total = aa_fluxes.sum(axis=0)







			for idx, massType in enumerate(massNames):
				massToPlot = mass.readColumn(massNames[idx])
				axesList[idx].plot(time / 60. / 60., massToPlot, linewidth = 2)

				axesList[idx].set_ylabel(cleanNames[idx] + " (fg)")

		for axes in axesList:
			axes.get_ylim()
			axes.set_yticks(list(axes.get_ylim()))

		axesList[0].set_title("Cell mass fractions")
		axesList[len(massNames) - 1].set_xlabel("Time (hr)")

		plt.subplots_adjust(hspace = 0.2, wspace = 0.5)
		fig.tight_layout()

		exportFigure(plt, plotOutDir, plotOutFileName,metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
