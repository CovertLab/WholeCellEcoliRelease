from __future__ import absolute_import, division, print_function

# constants for use as events in message passing

# environment and simulation control messages
TRIGGER_AGENT = 'TRIGGER_AGENT'
PAUSE_AGENT = 'PAUSE_AGENT'
SHUTDOWN_AGENT = 'SHUTDOWN_AGENT'
DIVIDE_CELL = 'DIVIDE_CELL'

# events for inner and outer communication
CELL_DECLARE = 'CELL_DECLARE'
CELL_INITIALIZE = 'CELL_INITIALIZE'
CELL_EXCHANGE = 'CELL_EXCHANGE'
CELL_SHUTDOWN = 'CELL_SHUTDOWN'
ENVIRONMENT_SYNCHRONIZE = 'ENVIRONMENT_SYNCHRONIZE'
ENVIRONMENT_UPDATE = 'ENVIRONMENT_UPDATE'

# events for the agent shepherd
ADD_AGENT = 'ADD_AGENT'
REMOVE_AGENT = 'REMOVE_AGENT'
TRIGGER_ALL = 'TRIGGER_ALL'
PAUSE_ALL = 'PAUSE_ALL'
SHUTDOWN_ALL = 'SHUTDOWN_ALL'

# universal agent messages
GLOBAL_SHUTDOWN = 'GLOBAL_SHUTDOWN'
