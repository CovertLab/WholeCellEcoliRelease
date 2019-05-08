'''
constants.py

Simulation constants.  Biological constants should go into the knowledge base;
configurable constants belong in the default_config.cfg file.
'''

from __future__ import absolute_import, division, print_function


SERIALIZED_RAW_DATA = "rawData.cPickle"
SERIALIZED_RAW_VALIDATION_DATA = "rawValidationData.cPickle"
SERIALIZED_VALIDATION_DATA = "validationData.cPickle"
SERIALIZED_SIM_DATA_FILENAME = "simData.cPickle"
SERIALIZED_SIM_DATA_MODIFIED = "simData_Modified.cPickle"
SERIALIZED_INHERITED_STATE = "Daughter%d_inherited_state.cPickle"

JSON_METADATA_FILE = 'metadata.json'

REQUEST_PRIORITY_DEGRADATION = 10
REQUEST_PRIORITY_DEFAULT = 0
REQUEST_PRIORITY_INTERN = -1 # processes that just request molecules
REQUEST_PRIORITY_METABOLISM = -10
REQUEST_PRIORITY_TWO_COMPONENT_SYSTEM = -5
