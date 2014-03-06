#!/usr/bin/env python

'''
config_test.py

Stores configuration information for tests
'''

import os
import json

# Number of processes to use in parallel tests
N_PROCESS = 4

with open(os.path.join('wholecell','utils','configfile','config_file_test.config')) as config_file:
	configs = json.loads(config_file.readline())
	N_PROCESSES = configs['N_PROCESSES']