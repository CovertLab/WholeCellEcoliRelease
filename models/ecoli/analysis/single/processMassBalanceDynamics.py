from __future__ import absolute_import, division, print_function

import cPickle
from itertools import izip
import os

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

THRESHOLD = 1e-13  # roughly, the mass of an electron


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)

		main_reader = TableReader(os.path.join(simOutDir, "Main"))
		mass = TableReader(os.path.join(simOutDir, "Mass"))
		fba_results = TableReader(os.path.join(simOutDir, "FBAResults"))

		initialTime = main_reader.readAttribute("initialTime")
		time = (main_reader.readColumn("time") - initialTime) / 60.
		n_timesteps = time.shape[0]

		processMassDifferences = mass.readColumn("processMassDifferences")
		processNames = mass.readAttribute("processNames")
		metabolism_index = processNames.index('Metabolism')

		# Adjust metabolism for exchange fluxes
		## Exchange fluxes will not capture rounding to single molecule level
		## so need to use change in metabolites as the mass coming into the cell
		metabolites = fba_results.readAttribute('outputMoleculeIDs')
		delta_metabolites = fba_results.readColumn('deltaMetabolites')

		conversion = 1e15 / sim_data.constants.nAvogadro.asNumber(1 / units.mol)
		metabolism_mass_imported = np.zeros(delta_metabolites.shape[0])
		for mol, flux in izip(metabolites, delta_metabolites.T):
			mol_mass = sim_data.getter.getMass([mol])[0].asNumber(units.g / units.mol)
			metabolism_mass_imported += mol_mass * conversion * flux

		## Normalize with initial mass to get relative mass difference
		metabolism_mass_difference = processMassDifferences[:, metabolism_index]
		adjusted_metabolism = (metabolism_mass_difference - metabolism_mass_imported).reshape(-1, 1)
		processMassDifferences = np.hstack((
			processMassDifferences[:, :metabolism_index+1].reshape(n_timesteps, -1),
			adjusted_metabolism,
			processMassDifferences[:, metabolism_index+1:].reshape(n_timesteps, -1),
			))
		processNames =  processNames[:metabolism_index+1] + ['Metabolism\n(exchange adjusted)'] + processNames[metabolism_index+1:]
		processMassDifferences[processMassDifferences == 0] = np.nan

		plt.figure(figsize = (8.5, 11))

		n_processes = len(processNames)

		n_cols = int(np.sqrt(n_processes))
		n_rows = int(np.ceil(n_processes/n_cols))

		axis = [time.min(), time.max(), np.nanmin(np.abs(processMassDifferences)), np.nanmax(np.abs(processMassDifferences))]

		for i, processName in enumerate(processNames):
			plt.subplot(n_rows, n_cols, i+1)

			series = np.abs(processMassDifferences[:, i])
			t = time.copy()

			t = t[series != 0]
			series = series[series != 0]

			if np.any(series > 0):
				colors = np.zeros((series.size, 3), np.float64)

				colors[np.abs(series) > THRESHOLD, :] = (1.0, 0.0, 0.0)

				# markers = ["^" if value > 0 else "v" for value in series]

				plt.scatter(
					t,
					np.abs(series),
					color = colors,
					marker = '.'
					# marker = markers
					)

			plt.hlines(THRESHOLD,  axis[0], axis[1], "r", "dashed")

			plt.title(processName)

			# plt.xlabel("time (min)")
			# plt.ylabel("mass diff ($|(m_f - m_i)/(m_i)|$)")

			plt.yscale("log")

			plt.axis(axis)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
