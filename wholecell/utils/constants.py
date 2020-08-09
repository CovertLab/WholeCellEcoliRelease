'''
constants.py

Simulation constants.  Biological constants should go into the knowledge base;
configurable constants belong in the default_config.cfg file.

'''

import os

TEST_FIXTURE_DIR = os.path.join("fixtures", "test")

OUTPUT_DIRECTORY = os.path.join("out", "simOut")

SERIALIZED_KB_DIR = os.path.join("fixtures", "kb")
SERIALIZED_KB_PREFIX = "KnowledgeBase"
SERIALIZED_KB_SUFFIX = ".pickle"
SERIALIZED_KB_MOST_FIT_FILENAME = SERIALIZED_KB_PREFIX + "_Most_Fit" + SERIALIZED_KB_SUFFIX
SERIALIZED_KB_UNFIT_FILENAME = SERIALIZED_KB_PREFIX + "_Unfit" + SERIALIZED_KB_SUFFIX

SERIALIZED_RAW_DATA = "rawData.pickle"
SERIALIZED_RAW_VALIDATION_DATA = "rawValidationData.pickle"
SERIALIZED_VALIDATION_DATA = "validationData.pickle"
SERIALIZED_SIM_DATA_PREFIX = "simData"
SERIALIZED_SIM_DATA_SUFFIX = ".pickle"
SERIALIZED_SIM_DATA_FILENAME = "simData_Fit_1.pickle"
SERIALIZED_FIT1_FILENAME = SERIALIZED_SIM_DATA_PREFIX + "_Fit_1" + SERIALIZED_SIM_DATA_SUFFIX
SERIALIZED_SIM_DATA_MOST_FIT_FILENAME = SERIALIZED_SIM_DATA_PREFIX + "_Most_Fit" + SERIALIZED_SIM_DATA_SUFFIX
SERIALIZED_SIM_DATA_MODIFIED = "simData_Modified.pickle"
SERIALIZED_METADATA_FILE = 'metadata.pickle'
SERIALIZED_CELL_SPECS = "cellSpecs.pickle"

REQUEST_PRIORITY_DEGRADATION = 10
REQUEST_PRIORITY_ATP_USAGE = 1
REQUEST_PRIORITY_DEFAULT = 0
REQUEST_PRIORITY_INTERN = -1 # processes that just request molecules
REQUEST_PRIORITY_METABOLISM = -10
REQUEST_PRIORITY_TWO_COMPONENT_SYSTEM = -5
