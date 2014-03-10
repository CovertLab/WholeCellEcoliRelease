#!/usr/bin/env python

'''
config.py

Stores configuration information
'''

import os
import json
from numpy import array,matrix

with open(os.path.join('wholecell','utils','configfile','config_file.config')) as config_file:
	configs = json.loads(config_file.readline())
	# Takes in keys and turns them into variables with same values as in dict:
	for k in configs:
		exec('{KEY} = {VALUE}'.format(KEY = k, VALUE = repr(configs[k])))