"""
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/27/2014
"""

from __future__ import absolute_import, division, print_function


import os

import numpy as np
from matplotlib import pyplot as plt
import six
from six.moves import cPickle, zip

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import units


THRESHOLD = 1e-13 # roughly, the mass of an electron

FG_PER_DALTON = 1.6605402e-9

# TODO: get these from the KB
REPRESENTATIVE_MASSES = {
	"proton":1.007 * FG_PER_DALTON,
	"amino acid":109 * FG_PER_DALTON,
	"ATP":551 * FG_PER_DALTON,
	"protein":40e3 * FG_PER_DALTON,
	"ribosome":2700e3 * FG_PER_DALTON
	}


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)

		# Listeners used
		main_reader = TableReader(os.path.join(simOutDir, "Main"))
		mass = TableReader(os.path.join(simOutDir, "Mass"))
		fba_results = TableReader(os.path.join(simOutDir, "FBAResults"))

		# Sim time
		initialTime = main_reader.readAttribute("initialTime")
		time = main_reader.readColumn("time") - initialTime

		# Mass differences
		processNames = mass.readAttribute("processNames")
		processMassDifferences = mass.readColumn("processMassDifferences")

		# Adjust metabolism for exchange fluxes
		## Exchange fluxes will not capture rounding to single molecule level
		## so need to use change in metabolites as the mass coming into the cell
		metabolites = fba_results.readAttribute('outputMoleculeIDs')
		delta_metabolites = fba_results.readColumn('deltaMetabolites')

		conversion = 1e15 / sim_data.constants.n_Avogadro.asNumber(1 / units.mol)
		metabolism_mass_imported = np.zeros(delta_metabolites.shape[0])
		for mol, flux in zip(metabolites, delta_metabolites.T):
			mol_mass = sim_data.getter.get_mass([mol])[0].asNumber(units.g / units.mol)
			metabolism_mass_imported += mol_mass * conversion * flux

		metabolism_mass_difference = processMassDifferences[:, processNames.index('Metabolism')]
		adjusted_metabolism = (metabolism_mass_difference - metabolism_mass_imported).reshape(-1, 1)
		processMassDifferences = np.hstack((adjusted_metabolism, processMassDifferences))

		# Average differences over cell cycle
		processNames = ['Metabolism\n(exchange adjusted)'] + processNames
		avgProcessMassDifferences = np.abs(processMassDifferences).sum(axis = 0) / len(time)
		index = np.arange(len(processNames))
		width = 1

		plt.figure(figsize = (8.5, 11))

		axes = plt.axes()

		r1 = axes.barh(index, avgProcessMassDifferences * (avgProcessMassDifferences > THRESHOLD), width, log = True, color = (0.9, 0.2, 0.2))
		r2 = axes.barh(index, avgProcessMassDifferences * (avgProcessMassDifferences <= THRESHOLD), width, log = True, color = (0.2, 0.2, 0.9))

		axes.set_yticks(index)
		axes.set_yticks(index + width/2, minor=True)
		axes.set_yticklabels(processNames) #, rotation = -45)

		axes.plot([THRESHOLD, THRESHOLD], [index[0], index[-1]+width], 'k--', linewidth=3)

		plt.text(THRESHOLD, index[-1], "electron", rotation = "vertical", va = "center", ha = "right")

		for name, mass in six.viewitems(REPRESENTATIVE_MASSES):
			plt.axvline(mass, color = "k")
			plt.text(mass, index[-1], name, rotation = "vertical", va = "center", ha = "right")

		plt.xlabel("Mass difference (fg)")

		plt.title("Average absolute change in mass by individual processes")

		plt.tight_layout()
		plt.grid(True, which="major", axis='x')
		plt.grid(True, which="minor", axis='y')

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
