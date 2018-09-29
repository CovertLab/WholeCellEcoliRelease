from __future__ import absolute_import, division, print_function

# constants for use as events in message passing

# environment and simulation control messages
TRIGGER_AGENT = 'TRIGGER_AGENT'
PAUSE_AGENT = 'PAUSE_AGENT'
DIVIDE_CELL = 'DIVIDE_CELL'
GLOBAL_SHUTDOWN = 'GLOBAL_SHUTDOWN'

# events for inner and outer communication
CELL_INITIALIZE = 'CELL_INITIALIZE'
CELL_EXCHANGE = 'CELL_EXCHANGE'
CELL_DIVISION = 'CELL_DIVISION'
CELL_SHUTDOWN = 'CELL_SHUTDOWN'
ENVIRONMENT_SYNCHRONIZE = 'ENVIRONMENT_SYNCHRONIZE'
ENVIRONMENT_UPDATE = 'ENVIRONMENT_UPDATE'

# events for the agent shepherd
ADD_AGENT = 'ADD_AGENT'
REMOVE_AGENT = 'REMOVE_AGENT'
SHUTDOWN_AGENT = 'SHUTDOWN_AGENT'
