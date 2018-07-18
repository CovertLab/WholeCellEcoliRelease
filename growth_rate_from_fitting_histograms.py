
# TODO (John): Move this script to a more appropriate place.

from __future__ import division

import os.path as pa

import numpy as np

import matplotlib.pyplot as plt

from wholecell.io.tablereader import TableReader

# Data sources
# TODO (John): Acquire these in a more programmatic fashion?
# TODO (John): os.path.join
DIRECTORY_NO_RIB_NO_RNAP = 'out/20180717.113833.402193__Neither_fit'
DIRECTORY_YES_RIB_NO_RNAP = None
DIRECTORY_NO_RIB_YES_RNAP = None
DIRECTORY_YES_RIB_YES_RNAP = None

DIRECTORIES = (
	DIRECTORY_NO_RIB_NO_RNAP,
	DIRECTORY_YES_RIB_NO_RNAP,
	DIRECTORY_NO_RIB_YES_RNAP,
	DIRECTORY_YES_RIB_YES_RNAP
	)

TITLES = ( # TODO (John): check order
	'({}) ribosomes, ({}) RNA polymerases'.format(rib, rnap)
	for rib in ('-', '+')
	for rnap in ('-', '+')
	)

N_SEEDS = 1 # TODO (John): 8
USE_GEN = (True, True, True, True) # TODO (John): only retain last two
N_GENS = len(USE_GEN)

# Plotting details
FIGSIZE = (8, 8) # size in inches
TARGET_DOUBLING_TIME_MINUTES = 45 # TODO (John): obtain from sim_data?
N_BINS = 30

TICK_FONTSIZE = 10
TITLE_FONTSIZE = 12

HIST_STYLE = dict(
	fc = 'royalblue'
	)

TARGET_LINE_STYLE = dict(
	color = 'crimson',
	lw = 2
	)

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

			tr = TableReader(path_to_table)

			time = tr.readColumn('time')

			doubling_time = (time[-1] - time[0]) / 60

			doubling_times.append(doubling_time)

	return doubling_times

plt.figure(figsize = FIGSIZE)

doubling_times = []

for directory in DIRECTORIES:
	if directory is None:
		doubling_times.append([])

	else:
		doubling_times.append(load_doubling_times(directory))

all_doubling_times = np.concatenate(doubling_times)

MIN_YMAX = 1
YMAX_FRAC_OVERSHOOT = 0.1
MIN_FRAC_SPREAD = 0.1

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

for (subplot_index, (dts, title)) in enumerate(zip(doubling_times, TITLES)):
	plt.subplot(2, 2, subplot_index + 1)

	freq = plt.hist(dts, bins = bins, **HIST_STYLE)[0]

	plt.xlim(bins[0], bins[-1])
	plt.ylim(0, max(np.max(freq), MIN_YMAX) * (1 + YMAX_FRAC_OVERSHOOT))

	plt.xticks(fontsize = TICK_FONTSIZE)
	plt.yticks(fontsize = TICK_FONTSIZE)

	plt.axvline(TARGET_DOUBLING_TIME_MINUTES, **TARGET_LINE_STYLE)

	plt.title(title, fontsize = TITLE_FONTSIZE)

plt.tight_layout()

plt.savefig('temp.pdf') # TODO (John): proper output path

# TODO (John): if __name__ == '__main__'
