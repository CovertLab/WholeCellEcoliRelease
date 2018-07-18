
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

N_SEEDS = 8
USE_GEN = [True, True, True, True] # TODO (John): only retain last two
N_GENS = len(USE_GEN)

# Plotting details
FIGSIZE = (8, 8) # size in inches
TARGET_DOUBLING_TIME_MINUTES = 45 # TODO (John): obtain from sim_data?

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

			time = tr.readColumn('Time')

			doubling_time = time[-1] - time[0]

			doubling_times.append(doubling_time)

	return doubling_times

def plot_doubling_time_distribution(doubling_times):
	pass

out = load_doubling_times(DIRECTORY_NO_RIB_NO_RNAP)

import ipdb; ipdb.set_trace()

# TODO (John): if __name__ == '__main__'
