'''
constants.py

Simulation constants.  Biological constants should go into the knowledge base;
configurable constants belong in the default_config.cfg file.

'''

import os

TEST_FIXTURE_DIR = os.path.join("fixtures", "test")

OUTPUT_DIRECTORY = os.path.join("out", "simOut")

SERIALIZED_KB_DIR = os.path.join("fixtures", "kb")
SERIALIZED_KB_FIT_FILENAME = "KnowledgeBase_Fit.cPickle"
SERIALIZED_KB_UNFIT_FILENAME = "KnowledgeBase_Unfit.cPickle"

REQUEST_PRIORITY_DEGRADATION = 10
REQUEST_PRIORITY_DEFAULT = 0
REQUEST_PRIORITY_INTERN = -1 # processes that just request molecules
REQUEST_PRIORITY_METABOLISM = -10
