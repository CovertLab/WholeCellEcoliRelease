"""
Plots for interpolation functions.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/8/20
"""

from __future__ import absolute_import, division, print_function

import os

from matplotlib import pyplot as plt
from matplotlib import gridspec
from six.moves import cPickle
import numpy as np

from models.ecoli.analysis import parcaAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.utils import constants, units


def get_raw(data, x_col, y_col, factor=1):
	xs = []
	ys = []
	for row in data:
		xs.append(row[x_col].asNumber(units.min))
		y = row[y_col] * factor
		if units.hasUnit(y):
			y = y.asNumber()
		ys.append(y)
	return xs, ys


class Plot(parcaAnalysisPlot.ParcaAnalysisPlot):
	def do_plot(self, input_dir, plot_out_dir, plot_out_filename, sim_data_file, validation_data_file, metadata):
		with open(os.path.join(input_dir, constants.SERIALIZED_RAW_DATA), 'rb') as f:
			raw_data = cPickle.load(f)
		with open(sim_data_file, 'rb') as f:
			sim_data = cPickle.load(f)
		growth = sim_data.growthRateParameters
		mass = sim_data.mass

		# Mapping of functions that perform interpolation to raw data
		interpolation_functions = {
			(growth.getFractionActiveRibosome, None):
				get_raw(raw_data.growthRateDependentParameters, 'doublingTime',
					'fractionActiveRibosome'),
			(growth.getFractionActiveRnap, None):
				get_raw(raw_data.growthRateDependentParameters, 'doublingTime',
					'fractionActiveRnap'),
			(growth.getppGppConc, None):
				get_raw(raw_data.growthRateDependentParameters, 'doublingTime',
					'ppGpp_conc', factor=growth._per_dry_mass_to_per_volume),
			(growth.getRibosomeElongationRate, None):
				get_raw(raw_data.growthRateDependentParameters, 'doublingTime',
					'ribosomeElongationRate'),
			(growth.getRnapElongationRate, None):
				get_raw(raw_data.growthRateDependentParameters, 'doublingTime',
					'rnaPolymeraseElongationRate'),
			(mass.get_dna_critical_mass, None): None,
			(mass.getAvgCellDryMass, None):
				get_raw(raw_data.dryMassComposition, 'doublingTime',
					'averageDryMass'),
			}

		# Interpolation functions that return values in a dictionary
		interpolation_functions.update({
			(mass.getMassFraction, fraction):
				get_raw(raw_data.dryMassComposition, 'doublingTime',
					'{}MassFraction'.format(fraction))
			for fraction in mass.getMassFraction(45*units.min)
			})
		interpolation_functions.update({
			(mass.getFractionMass, fraction): None
			for fraction in mass.getFractionMass(45*units.min)
			})

		# TODO: handle getTrnaDistribution and all 86 outputs from 'molar_ratio_to_16SrRNA'
		# mass.getTrnaDistribution(45.*units.min)['molar_ratio_to_16SrRNA']

		# Doubling times to show on plot (extended range including all conditions)
		doubling_times = np.unique([
			dt.asNumber(units.min)
			for dt in sim_data.conditionToDoublingTime.values()
			])
		doubling_time_range = np.arange(0.5 * doubling_times.min(), 1.2 * doubling_times.max())

		# Create Plot
		plt.figure(figsize=(20, 20))
		n_plots = len(interpolation_functions)
		cols = 5
		gs = gridspec.GridSpec(int(np.ceil(n_plots / cols)), cols)

		for i, ((fun, key), data) in enumerate(interpolation_functions.items()):
			ax = plt.subplot(gs[i // cols, i % cols])

			# Get interpolation values and handle units
			values = []
			unit = None
			for dt in units.min * doubling_time_range:
				try:
					value = fun(dt)
					if key:
						value = value[key]

					if units.hasUnit(value):
						unit = str(units.getUnit(value)).strip('1 []')
						value = value.asNumber()
				except Exception as e:
					value = np.nan
				values.append(value)
			y_interp = np.array(values)

			# Create y labels
			label = fun.__name__ + ('\n{}'.format(key) if key else '')
			if unit:
				label += '\n({})'.format(unit)

			# Plot data
			ax.plot(doubling_time_range, y_interp)
			if data:
				ax.plot(data[0], data[1], 'or')
			for dt in doubling_times:
				ax.axvline(dt, linestyle='--', color='k', linewidth=0.5)

			# Formatting
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)
			ax.set_ylabel(label)
			ax.set_xlim([np.min(doubling_time_range), np.max(doubling_time_range)])

		## Save figure
		plt.tight_layout()
		exportFigure(plt, plot_out_dir, plot_out_filename, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
