#!/usr/bin/env python

'''
config_fixtures.py

Stores configuration information for simulation
'''

import os
import json

# Knowledge base direcotry (default configuration)
KNOWLEDGEBASE_DIRECTORY = '~/Documents/kbEcoli/'

with open(os.path.join('wholecell','utils','configfile','config_file_sim.config')) as config_file:
	configs = json.loads(config_file.readline())
	KNOWLEDGEBASE_DIRECTORY = configs['KNOWLEDGEBASE_DIRECTORY']