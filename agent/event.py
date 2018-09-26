from __future__ import absolute_import, division, print_function

# constants for use as events in message passing

# events for inner and outer communication
TRIGGER_EXECUTION = 'TRIGGER_EXECUTION'
SIMULATION_INITIALIZED = 'SIMULATION_INITIALIZED'
SIMULATION_ENVIRONMENT = 'SIMULATION_ENVIRONMENT'
SIMULATION_SHUTDOWN = 'SIMULATION_SHUTDOWN'
SYNCHRONIZE_SIMULATION = 'SYNCHRONIZE_SIMULATION'
ENVIRONMENT_UPDATED = 'ENVIRONMENT_UPDATED'
DIVIDE_CELL = 'DIVIDE_CELL'
CELL_DIVISION = 'CELL_DIVISION'
PAUSE_ENVIRONMENT = 'PAUSE_ENVIRONMENT'
GLOBAL_SHUTDOWN = 'GLOBAL_SHUTDOWN'

# events for the agent shepherd
ADD_AGENT = 'ADD_AGENT'
REMOVE_AGENT = 'REMOVE_AGENT'
SHUTDOWN_AGENT = 'SHUTDOWN_AGENT'
