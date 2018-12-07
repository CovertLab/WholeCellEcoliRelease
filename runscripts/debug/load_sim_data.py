'''
Script to quickly load raw_data and sim_data to interactively explore contents
through ipdb instead of running the fitter. Shortcut names rd and sd are created
for typing ease.

Use:
	python runscripts/manual/load_sim_data.py

Requires:
	cached/ directory in wcEcoli root with rawData.cPickle and simData_Fit_1.cPickle

Tip:
Create an alias in your .bash_profile to easily run this by just typing lsd (or your
alias of choice):
	alias lsd="python runscripts/manual/load_sim_data.py"
'''

import cPickle
import os

import ipdb
import numpy as np


root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
data_dir = os.path.join(root_dir, 'cached')
raw_data_filename = 'rawData.cPickle'
sim_data_filename = 'simData_Fit_1.cPickle'

raw_data = cPickle.load(open(os.path.join(data_dir, raw_data_filename), 'rb'))
sim_data = cPickle.load(open(os.path.join(data_dir, sim_data_filename), 'rb'))

rd = raw_data
sd = sim_data

ipdb.set_trace()
