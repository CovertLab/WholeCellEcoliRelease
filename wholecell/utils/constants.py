'''
constants.py

Simulation constants.  Biological constants should go into the knowledge base.
'''

from __future__ import absolute_import, division, print_function


OPERON_OPTIONS = ('off', 'on')  # TODO(jerry): add 'both'
DEFAULT_OPERON_OPTION = 'off'

SERIALIZED_RAW_DATA = "rawData.cPickle"
SERIALIZED_RAW_VALIDATION_DATA = "rawValidationData.cPickle"
SERIALIZED_VALIDATION_DATA = "validationData.cPickle"
SERIALIZED_SIM_DATA_FILENAME = "simData.cPickle"
SERIALIZED_METRICS_DATA_FILENAME = "metricsData.cPickle"
SERIALIZED_SIM_DATA_MODIFIED = "simData_Modified.cPickle"
SERIALIZED_INHERITED_STATE = "Daughter%d_inherited_state.cPickle"

# Workflow directories
# TODO: add 'plotOut', etc.
KB_PLOT_OUTPUT_DIR = 'kb_plot_out'
KB_DIR = 'kb'  # TODO: replace hardcoded 'kb' strings with this variable in codebase

JSON_METADATA_FILE = 'metadata.json'
GIT_DIFF_FILE = 'git_diff.txt'

REQUEST_PRIORITY_DEGRADATION = 10
REQUEST_PRIORITY_DEFAULT = 0
REQUEST_PRIORITY_INTERN = -1 # processes that just request molecules
REQUEST_PRIORITY_TWO_COMPONENT_SYSTEM = -5
REQUEST_PRIORITY_TF_BINDING = -10  # need to have low priority with requestAll
REQUEST_PRIORITY_METABOLISM = -10
