from __future__ import absolute_import, division, print_function

# constants for use as events in message passing

# environment and simulation control messages
TRIGGER_AGENT = 'TRIGGER_AGENT'  # Start/resume simulation.
PAUSE_AGENT = 'PAUSE_AGENT'
SHUTDOWN_AGENT = 'SHUTDOWN_AGENT'
DIVIDE_CELL = 'DIVIDE_CELL'  # Ask a cell to divide now.

# events from inner to outer
CELL_DECLARE = 'CELL_DECLARE'  # A new cell agent has an ID and can be displayed but its simulation is still initializing.
CELL_INITIALIZE = 'CELL_INITIALIZE'  # A new cell agent is initialized and ready to run.
CELL_EXCHANGE = 'CELL_EXCHANGE'  # Molecule exchange data and maybe cell division data.
CELL_SHUTDOWN = 'CELL_SHUTDOWN'  # A cell is shutting down.

# events from outer to inner
ENVIRONMENT_SYNCHRONIZE = 'ENVIRONMENT_SYNCHRONIZE'  # Synchronize the cell's clock to the environment.
ENVIRONMENT_UPDATE = 'ENVIRONMENT_UPDATE'  # Molecule exchange data.

# events to an agent shepherd
ADD_AGENT = 'ADD_AGENT'  # Spawn an agent subprocess.
REMOVE_AGENT = 'REMOVE_AGENT'
TRIGGER_ALL = 'TRIGGER_ALL'
PAUSE_ALL = 'PAUSE_ALL'
SHUTDOWN_ALL = 'SHUTDOWN_ALL'

# universal agent messages
GLOBAL_SHUTDOWN = 'GLOBAL_SHUTDOWN'  # Shut down all agents.
