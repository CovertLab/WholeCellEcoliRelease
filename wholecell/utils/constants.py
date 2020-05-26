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
SERIALIZED_METRICS_DATA_FILENAME = "metricsData.cPickle"
SERIALIZED_SIM_DATA_MODIFIED = "simData_Modified.cPickle"
SERIALIZED_INHERITED_STATE = "Daughter%d_inherited_state.cPickle"

# Workflow directories
# TODO: add 'plotOut', 'kb', etc.
KB_PLOT_OUTPUT_DIR = 'kb_plot_out'

JSON_METADATA_FILE = 'metadata.json'

REQUEST_PRIORITY_DEGRADATION = 10
REQUEST_PRIORITY_DEFAULT = 0
REQUEST_PRIORITY_INTERN = -1 # processes that just request molecules
REQUEST_PRIORITY_TWO_COMPONENT_SYSTEM = -5
REQUEST_PRIORITY_TF_BINDING = -10  # need to have low priority with requestAll
REQUEST_PRIORITY_METABOLISM = -10
