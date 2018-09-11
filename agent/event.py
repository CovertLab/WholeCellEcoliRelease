from __future__ import absolute_import, division, print_function

# constants for use as events in message passing

# events for inner and outer communication
TRIGGER_EXECUTION = 'TRIGGER_EXECUTION'
SIMULATION_INITIALIZED = 'SIMULATION_INITIALIZED'
SIMULATION_ENVIRONMENT = 'SIMULATION_ENVIRONMENT'
SIMULATION_SHUTDOWN = 'SIMULATION_SHUTDOWN'
SYNCHRONIZE_SIMULATION = 'SYNCHRONIZE_SIMULATION'
ENVIRONMENT_UPDATED = 'ENVIRONMENT_UPDATED'
SHUTDOWN_ENVIRONMENT = 'SHUTDOWN_ENVIRONMENT'
SHUTDOWN_SIMULATION = 'SHUTDOWN_SIMULATION'
GLOBAL_SHUTDOWN = 'GLOBAL_SHUTDOWN'

# events for the agent shepherd
INITIALIZE_AGENT = 'INITIALIZE_AGENT'
SHUTDOWN_AGENT = 'SHUTDOWN_AGENT'
