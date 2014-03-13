#!/usr/bin/env python

'''
config.py

Stores configuration information
'''

import os
import json
from numpy import array,matrix

hudson_config_file = os.path.join('wholecell','utils','configfile','config_file.config_hudson')
user_config_file = os.path.join('wholecell','utils','configfile','config_file.config')

if os.path.exists(user_config_file):
	with open(user_config_file) as config_file:
		configs = json.loads(config_file.readline())
		# Takes in keys and turns them into variables with same values as in dict:
		for k in configs:
			exec('{KEY} = {VALUE}'.format(KEY = k, VALUE = repr(configs[k])))
else:
	import warnings
	warnings.warn('Using Hudson config file! To create your own copy the config_file.config_default file.\n')
	
	with open(user_config_file) as config_file:
		configs = json.loads(config_file.readline())
		# Takes in keys and turns them into variables with same values as in dict:
		for k in configs:
			exec('{KEY} = {VALUE}'.format(KEY = k, VALUE = repr(configs[k])))