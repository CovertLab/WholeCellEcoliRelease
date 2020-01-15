"""
Plots rnap capacity

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/18/2015
"""

from __future__ import division

import cPickle

from matplotlib import pyplot as plt
import numpy as np
import os

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.analysis.analysis_tools import read_bulk_molecule_counts
from wholecell.io.tablereader import TableReader
from wholecell.utils import filepath
from wholecell.utils import units

FONT = {
	'size':	8
	}


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		filepath.makedirs(plotOutDir)

		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)

		# Load data from KB
		nAvogadro = sim_data.constants.nAvogadro

		# Listeners used
		unique_molecules_reader = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))
		main_reader = TableReader(os.path.join(simOutDir, "Main"))
		rnap_reader = TableReader(os.path.join(simOutDir, "RnapData"))

		# Get ID and mass of inactive (bulk) RNAP
		inactive_rnap_id = sim_data.moleculeIds.rnapFull
		rnap_mass = sim_data.getter.getMass([inactive_rnap_id])[0]

		# Read time data
		initial_time = main_reader.readAttribute("initialTime")
		time = main_reader.readColumn("time") - initial_time
		timeStep = main_reader.readColumn("timeStepSec")

		# Calculate the elongation rate for the given condition
		nutrients = sim_data.conditions[sim_data.condition]["nutrients"]
		elongation_rate = sim_data.process.transcription.rnaPolymeraseElongationRateDict[nutrients].asNumber(units.aa/units.s) * timeStep

		# Load rnap data
		actual_elongations = rnap_reader.readColumn("actualElongations")

		# Load counts of subunits and active ribosomes
		(inactive_rnap_counts, ) = read_bulk_molecule_counts(
			simOutDir, ([inactive_rnap_id], ))
		active_rnap_index = unique_molecules_reader.readAttribute("uniqueMoleculeIds").index('active_RNAP')
		active_rnap_counts = unique_molecules_reader.readColumn("uniqueMoleculeCounts")[:, active_rnap_index]

		# Calculate statistics
		total_rnap_counts = active_rnap_counts + inactive_rnap_counts
		total_rnap_capacity = total_rnap_counts * elongation_rate

		mass_fraction_active = active_rnap_counts/total_rnap_counts

		plt.figure(figsize = (8.5, 15))
		plt.rc('font', **FONT)

		total_rnap_capacity_axis = plt.subplot(6,1,1)
		total_rnap_capacity_axis.plot(
			time / 60., total_rnap_capacity,
			label="Theoretical total RNAP capacity", linewidth=2, color='b')
		total_rnap_capacity_axis.plot(
			time / 60., actual_elongations,
			label="Actual elongations", linewidth=2, color='r')
		total_rnap_capacity_axis.set_ylabel("Nucleotides polymerized")
		total_rnap_capacity_axis.legend(ncol=2)

		active_rnap_capacity_axis = plt.subplot(6,1,2)
		active_rnap_capacity_axis.plot(
			time / 60., active_rnap_counts * elongation_rate,
			label="Theoretical active RNAP capacity", linewidth=2, color='b')
		active_rnap_capacity_axis.plot(
			time / 60., actual_elongations,
			label="Actual elongations", linewidth=2, color='r')
		active_rnap_capacity_axis.set_ylabel("Nucleotides polymerized")
		active_rnap_capacity_axis.legend(ncol=2)

		inactive_rnap_capacity_axis = plt.subplot(6,1,3)
		inactive_rnap_capacity_axis.plot(
			time / 60., inactive_rnap_counts * elongation_rate,
			label="Theoretical inactive ribosome capacity", linewidth=2, color='b')
		inactive_rnap_capacity_axis.set_ylabel("Nucleotides polymerized")
		inactive_rnap_capacity_axis.legend(ncol=2)

		fractionalCapacity_axis = plt.subplot(6,1,4)
		fractionalCapacity_axis.plot(
			time / 60., actual_elongations / total_rnap_capacity,
			linewidth=2, color='k')
		fractionalCapacity_axis.set_ylabel("Fraction of total RNAP capacity used")

		effectiveElongationRate_axis = plt.subplot(6,1,5)
		effectiveElongationRate_axis.plot(
			time / 60., actual_elongations / active_rnap_counts,
			linewidth=2, color='k')
		effectiveElongationRate_axis.set_ylabel("Effective elongation rate (aa/s/rnap)")

		fractionActive_axis = plt.subplot(6,1,6)
		fractionActive_axis.plot(
			time / 60., mass_fraction_active,
			linewidth=2, color='k')
		fractionActive_axis.set_ylabel("Mass fraction of active RNAPs")
		fractionActive_axis.set_yticks(np.arange(0., 1.1, 0.1))

		# Save
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
