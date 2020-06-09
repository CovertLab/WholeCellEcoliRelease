'''
Script to quickly load raw_data and sim_data to interactively explore contents
through ipdb instead of running the parameter calculator (parca).
Shortcut names rd and sd are created for typing ease as well as variables for
each process object.

Use:
	python runscripts/debug/load_sim_data.py [path to simulation output root]

Requires:
	cached/ directory in wcEcoli root with rawData.cPickle and simData.cPickle
	if a simulation directory is not supplied

Tip:
Create an alias in your .bash_profile to easily run this by just typing lsd (or your
alias of choice):
	alias lsd="python runscripts/debug/load_sim_data.py"
'''

from __future__ import absolute_import, division, print_function

from six.moves import cPickle
import os
import sys

import ipdb
import numpy as np

from wholecell.utils import constants, filepath, units


# Imported just to make available in debugging; suppress the unused warning
_unused = [np, units]

# Parse the first arg as simulation root directory, use cached sim_data as default
if len(sys.argv) > 1:
	data_dir = os.path.join(sys.argv[1], 'kb')
else:
	data_dir = os.path.join(filepath.ROOT_PATH, 'cached')
raw_data_filename = constants.SERIALIZED_RAW_DATA
sim_data_filename = constants.SERIALIZED_SIM_DATA_FILENAME

raw_data = cPickle.load(open(os.path.join(data_dir, raw_data_filename), 'rb'))
sim_data = cPickle.load(open(os.path.join(data_dir, sim_data_filename), 'rb'))
print('Loaded {} from {}'.format(sim_data_filename, data_dir))

# Shortcuts for easier access or tab complete
rd = raw_data
sd = sim_data
process = sd.process
complexation = process.complexation
equilibrium = process.equilibrium
metabolism = process.metabolism
replication = process.replication
rna_decay = process.rna_decay
transcription = process.transcription
transcription_regulation = process.transcription_regulation
translation = process.translation
two_component_system = process.two_component_system

ipdb.set_trace()
