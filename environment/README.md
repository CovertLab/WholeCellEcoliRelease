# Environment

Models of environments, which can be interfaced with cells via the agent module. Includes boot.py for starting both
environments and cells.

## Setup

See the top-level [README.md](../README.md) for general setup instructions, and the
[agent README.md](../agent/README.md) for multi-agent simulation setup.

## Usage

1. To run a simulation, you will need to have one sim_data object for the cell. You can generate this object with the
runFitter manual runscript. In the wcEcoli directory:

    `> PYTHONPATH="$PWD" python runscripts/manual/runFitter.py`

2. You will also need to runSim in order to generate the metadata folder.

    `> PYTHONPATH="$PWD" python runscripts/manual/runSim.py`

3. With those files in place, you are ready to begin running the environmental simulation processes.
See [agent/README.md](../agent/README.md) for instructions to set up your Zookeeper and Kafka servers.

   1. Start Zookeeper in the directory where you untarred the Kafka and Zookeeper servers:

      `> ./bin/zookeeper-server-start.sh ./config/zookeeper.properties`

   2. Then start the Kafka server in another shell tab in the same directory:

      `> ./bin/kafka-server-start.sh ./config/server.properties`

4. Boot an environment agent and cell agents through this environment directory.
Boot each process in a new tab. (**Tip:** Use iTerm split windows to make
it easy to watch all these processeses at once.)

   1. In the first tab start an environment model:

      `> python -m environment.boot lattice`

      This starts the environment and opens a "microscope" view onto the plate, waiting for cell simulations to register.

   2. Now start a cell agent in a new tab:

      `> python -m environment.boot ecoli --id 1`

5. You can start as many cells as desired, each with its own unique id (agent name).
You will see a message sent from the newly initialized simulation on the `environment_listen` topic:

   `<-- environment_listen: {'event': 'SIMULATION_INITIALIZED', 'id': '1'}`

   and in the environment tab you will see it has received a message from the new cell simulation:

   `--> environment_listen: {u'event': u'SIMULATION_INITIALIZED', u'id': u'1'}`

6. Finally, use the agent module to trigger the environmental simulation in a separate "command" tab:

   `> python -m agent.boot trigger`

7. To stop the simulation, run `shutdown` in the command tab:

   `> python -m agent.boot shutdown`
