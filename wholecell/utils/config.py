#!/usr/bin/env python

'''
config.py

Stores configuration information
'''

from __future__ import absolute_import, division, print_function

import os
import json
import six

defaultConfig = json.load(
	open(os.path.join('wholecell', 'utils', 'default_config.cfg'))
	)

globals().update(defaultConfig)

# Load the user config file, or, if there is none, create it

userConfigPath = os.path.join('user', 'user_config.cfg')

if not os.path.exists(userConfigPath):
	print("Creating blank user configuration file at {}".format(userConfigPath))
	with open(userConfigPath, 'w') as configFile:
		json.dump({}, configFile)

userConfig = json.load(open(userConfigPath))

unknownOptions = six.viewkeys(userConfig) - six.viewkeys(defaultConfig)

if unknownOptions:
	raise Exception("Unknown configuration options defined in {}: {}".format(
		userConfigPath,
		', '.join(unknownOptions)
		))

globals().update(userConfig)
