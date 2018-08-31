# Environment

Models of environments, which can be interfaced with cells via the agent module. Includes boot.py for starting both
environments and cells.

## Usage
To run a simulation, you will need to have one sim_data object for the cell. You can generate this object with the
runFitter manual runscript. You will also need to runSim in order to generate the metadata folder.

With those files in place, you are ready to begin an environmental simulation. First, initiate the agent module for message
passing between the environment and cells. View agent/README.md for instructions of agent's requirements.
In brief, this includes:

Start Zookeeper in the zookeeper directory:

    > ./bin/zookeeper-server-start.sh ./config/zookeeper.properties

then start Kafka in the same directory:

    > ./bin/kafka-server-start.sh ./config/server.properties

Next, you will boot an environment and cells through this environment directory. Boot each process in a new tab.

In the first tab start an environment model:

    > python -m environment.boot lattice

This starts the environment, and is waiting for simulations initialize and register. Let's do that now. In a new tab:

    > python -m environment.boot ecoli --id 1

You can start as many cells as desired, each with its own unique id. Assigning a unique id to each simulation allows the
environment to communicate in specific ways with each simulation. You will see a message sent from the newly initialized
simulation on the `environment_listen` topic:

    <-- environment_listen: {'event': 'SIMULATION_INITIALIZED', 'id': '1'}

and also if you go back to the environment tab you will see it has received a message from the new simulation:

    --> environment_listen: {u'event': u'SIMULATION_INITIALIZED', u'id': u'1'}

Finally, use the agent module to trigger the environmental simulation in a separate `command` tab:

    > python -m agent.boot trigger

Once you are ready to stop the simulation, you can call shutdown from the command tab:

    > python -m agent.boot shutdown
