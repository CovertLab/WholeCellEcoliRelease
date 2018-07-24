
'''

This script generate a 2x2 figure of histograms, showing the doubling time for
simulations under alternative (limited) fitting procedures.  It can't be
implemetned as a standard analysis script, because the underlying feature
(passing options to the fitter) is not currently supported by our 'variant'
framework.

Usage
-----

This script should be executed from the root wcEcoli directory.  You will need
to adjust the directories for simulation output.

'''

from __future__ import division

import os.path as pa

import numpy as np

import matplotlib.pyplot as plt

from wholecell.io.tablereader import TableReader

# Data sources
# TODO (John): Acquire these in a more programmatic fashion?
# TODO (John): os.path.join
DIRECTORY_NO_RIB_NO_RNAP = 'out/20180717.160035.945718__Neither_fit'
DIRECTORY_YES_RIB_NO_RNAP = 'out/20180717.160011.414995__Ribosomes_fit'
DIRECTORY_NO_RIB_YES_RNAP = 'out/20180717.160023.903645__RNApoly_fit'
DIRECTORY_YES_RIB_YES_RNAP = 'out/20180717.155958.639346__Both_fit'

DIRECTORIES = (
	DIRECTORY_NO_RIB_NO_RNAP,
	DIRECTORY_YES_RIB_NO_RNAP,
	DIRECTORY_NO_RIB_YES_RNAP,
	DIRECTORY_YES_RIB_YES_RNAP
	)

TITLES = (
	'({}) ribosomes, ({}) RNA polymerases'.format(rib, rnap)
	for rnap in ('-', '+')
	for rib in ('-', '+')
	)

N_SEEDS = 8 # TODO (John): 8
USE_GEN = (True, True, True, True) # TODO (John): only retain last two
N_GENS = len(USE_GEN)

BOOTSTRAP = False # if True, adds more observations by distribution sampling
FOLD_BOOTSTRAP = 15 # number of additional observations to include

# Plotting details
FIGSIZE = (8, 8) # size in inches
TARGET_DOUBLING_TIME_MINUTES = 45 # TODO (John): obtain from sim_data?
N_BINS = 30
SHARE_YMAX = True

TICK_FONTSIZE = 10
TITLE_FONTSIZE = 12

HIST_STYLE = dict(
	fc = 'royalblue'
	)

TARGET_LINE_STYLE = dict(
	color = 'crimson',
	lw = 2
	)

MIN_YMAX = 1
YMAX_FRAC_OVERSHOOT = 0.1
MIN_FRAC_SPREAD = 0.1

def load_doubling_times(root_directory):
	doubling_times = []

	for seed in xrange(N_SEEDS):
		for (gen, do_use) in enumerate(USE_GEN):
			if not do_use: continue

			path_to_table = pa.join(
				root_directory,
				'wildtype_000000',
				'{:06n}'.format(seed),
				'generation_{:06n}'.format(gen),
				'0'*6,
				'simOut',
				'Main'
				)

			try:
				tr = TableReader(path_to_table)

				time = tr.readColumn('time')

			except Exception as e:
				print 'Failed to load {}'.format(path_to_table)
				continue

			doubling_time = (time[-1] - time[0]) / 60

			doubling_times.append(doubling_time)

	return doubling_times

def bootstrap(x, fold):
	n = x.size

	m = np.mean(x)
	s = np.std(x)

	new = m + s*np.random.normal(size = n*fold)

	return np.concatenate([x, new])

plt.figure(figsize = FIGSIZE)

doubling_times = []

for directory in DIRECTORIES:
	if directory is None:
		doubling_times.append([])

	else:
		dts = np.array(load_doubling_times(directory))

		if BOOTSTRAP:
			dts = bootstrap(dts, FOLD_BOOTSTRAP)

		doubling_times.append(dts)

all_doubling_times = np.concatenate(doubling_times)

smallest = min(
	all_doubling_times.min(),
	TARGET_DOUBLING_TIME_MINUTES * (1 - MIN_FRAC_SPREAD)
	)
largest = max(
	all_doubling_times.max(),
	TARGET_DOUBLING_TIME_MINUTES * (1 + MIN_FRAC_SPREAD)
	)

left = np.linspace(smallest, largest, N_BINS-1)
bins = np.concatenate([left, [left[1] - left[0] + left[-1]]])

max_freq = 0

for (subplot_index, (dts, title)) in enumerate(zip(doubling_times, TITLES)):
	plt.subplot(2, 2, subplot_index + 1)

	freq = plt.hist(dts, bins = bins, **HIST_STYLE)[0]

	max_freq = max(max_freq, freq.max())

	plt.xlim(bins[0], bins[-1])
	plt.ylim(0, max(np.max(freq), MIN_YMAX) * (1 + YMAX_FRAC_OVERSHOOT))

	plt.xticks(fontsize = TICK_FONTSIZE)
	plt.yticks(fontsize = TICK_FONTSIZE)

	plt.axvline(TARGET_DOUBLING_TIME_MINUTES, **TARGET_LINE_STYLE)

	plt.title(title + '\nn = {}'.format(dts.size), fontsize = TITLE_FONTSIZE)

if SHARE_YMAX:
	ymax = max(max_freq, MIN_YMAX) * (1 + YMAX_FRAC_OVERSHOOT)

	for subplot_index in xrange(len(doubling_times)):
		plt.subplot(2, 2, subplot_index + 1)
		plt.ylim(0, ymax)

plt.tight_layout()

plt.savefig('temp.png') # TODO (John): proper output path

# TODO (John): if __name__ == '__main__'
