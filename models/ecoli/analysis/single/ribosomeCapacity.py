"""
Plots ribosome capacity

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 11/20/2014
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
		ribosome_reader = TableReader(os.path.join(simOutDir, "RibosomeData"))

		# Get IDs of ribosome subunits
		ribosome_subunit_ids = [
			sim_data.moleculeIds.s50_fullComplex,
			sim_data.moleculeIds.s30_fullComplex,
			]

		# Get masses of full ribosomes and subunits
		ribosome_subunit_masses = sim_data.getter.getMass(ribosome_subunit_ids)
		full_ribosome_mass = units.sum(ribosome_subunit_masses)

		# Read time data
		initial_time = main_reader.readAttribute("initialTime")
		time = main_reader.readColumn("time") - initial_time
		timeStep = main_reader.readColumn("timeStepSec")

		# Calculate the elongation rate for the given condition
		nutrients = sim_data.conditions[sim_data.condition]["nutrients"]
		elongation_rate = sim_data.process.translation.ribosomeElongationRateDict[nutrients].asNumber(units.aa/units.s)

		# Load ribosome data
		actual_elongations = ribosome_reader.readColumn("actualElongations")
		actual_elongation_rate = actual_elongations / timeStep

		# Load counts of subunits and active ribosomes
		(ribosome_subunit_counts, ) = read_bulk_molecule_counts(
			simOutDir, (ribosome_subunit_ids, ))
		active_ribosome_index = unique_molecules_reader.readAttribute("uniqueMoleculeIds").index('active_ribosome')
		active_ribosome_counts = unique_molecules_reader.readColumn("uniqueMoleculeCounts")[:, active_ribosome_index]

		# Calculate statistics
		total_ribosome_counts = active_ribosome_counts + ribosome_subunit_counts.min(axis=1)
		total_ribosome_capacity = total_ribosome_counts * elongation_rate

		free_subunit_mass = (
			(ribosome_subunit_masses * ribosome_subunit_counts / nAvogadro).asNumber(units.fg)
			).sum(axis=1)
		active_ribosome_mass = (full_ribosome_mass * active_ribosome_counts / nAvogadro).asNumber(units.fg)
		total_ribosome_mass = free_subunit_mass + active_ribosome_mass
		mass_fraction_active = active_ribosome_mass / total_ribosome_mass

		plt.figure(figsize = (8.5, 15))
		plt.rc('font', **FONT)

		ribosomeCapacity_axis = plt.subplot(6,1,1)
		ribosomeCapacity_axis.plot(
			time / 60., total_ribosome_capacity,
			label="Theoretical total ribosome rate", linewidth=2, color='b')
		ribosomeCapacity_axis.plot(
			time / 60., actual_elongation_rate,
			label="Actual elongation rate", linewidth=2, color='r')
		ribosomeCapacity_axis.set_ylabel("Total amino acid\npolymerization rate\n(AA/s)")
		ribosomeCapacity_axis.legend(ncol=2)

		activeRibosomeCapacity_axis = plt.subplot(6,1,2)
		activeRibosomeCapacity_axis.plot(
			time / 60., active_ribosome_counts * elongation_rate,
			label="Theoretical active ribosome rate", linewidth=2, color='b')
		activeRibosomeCapacity_axis.plot(
			time / 60., actual_elongation_rate,
			label="Actual elongation rate", linewidth=2, color='r')
		activeRibosomeCapacity_axis.set_ylabel("Total amino acid\npolymerization rate\n(AA/s)")
		activeRibosomeCapacity_axis.legend(ncol=2)

		inactiveRibosomeCapacity_axis = plt.subplot(6,1,3)
		inactiveRibosomeCapacity_axis.plot(
			time / 60., ribosome_subunit_counts.min(axis=1) * elongation_rate,
			label="Theoretical inactive ribosome rate", linewidth=2, color='b')
		inactiveRibosomeCapacity_axis.set_ylabel("Total amino acid\npolymerization rate\n(AA/s)")
		inactiveRibosomeCapacity_axis.legend(ncol=2)

		fractionalCapacity_axis = plt.subplot(6,1,4)
		fractionalCapacity_axis.plot(
			time / 60., actual_elongation_rate / total_ribosome_capacity,
			linewidth=2, color='k')
		fractionalCapacity_axis.set_ylabel("Fraction of total ribosome capacity used")

		effectiveElongationRate_axis = plt.subplot(6,1,5)
		effectiveElongationRate_axis.plot(
			time / 60., actual_elongation_rate / active_ribosome_counts,
			linewidth=2, color='k')
		effectiveElongationRate_axis.set_ylabel("Relative elongation rate (aa/s/ribosome)")

		fractionActive_axis = plt.subplot(6,1,6)
		fractionActive_axis.plot(
			time / 60., mass_fraction_active,
			linewidth=2, color='k')
		fractionActive_axis.set_ylabel("Mass fraction of active ribosomes")
		fractionActive_axis.set_yticks(np.arange(0., 1.1, 0.1))

		# Save
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
