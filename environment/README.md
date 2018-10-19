# Environment

Models of environments, which can be interfaced with cells via the agent module. Includes boot.py for starting both
environments and cells.

## Setup

See the top-level [README.md](../README.md) for general setup instructions, and the
[agent README.md](../agent/README.md) for multi-agent simulation setup.

## Usage

1. To run Cell simulations, you need to have the sim_data files. You can generate them via the
runFitter manual runscript. In the wcEcoli directory:

    `> PYTHONPATH="$PWD" python runscripts/manual/runFitter.py`

2. See [agent/README.md](../agent/README.md) for instructions to set up your Zookeeper and Kafka servers. To recap:

   1. Start Zookeeper in the directory where you untarred the Kafka and Zookeeper software:

      `> ./bin/zookeeper-server-start.sh ./config/zookeeper.properties`

   2. Then start the Kafka server in another shell tab in the same directory:

      `> ./bin/kafka-server-start.sh ./config/server.properties`

3. **Optional:** Start the [Environment visualization](https://github.com/CovertLab/environment)
server and browser window per the instructions on that page. To recap:

   1. Run the visualization server in the root directory of that repository:

      `> lein run`

   2. Open a browser window onto [http://localhost:33332](http://localhost:33332)

4. Boot an Environment agent and Cell agents through this environment directory, with
each process in a new terminal tab. (**Tip:** Use iTerm split windows to make
it easy to watch all these processeses at once.)

   1. In the first tab start an Environment model:

      `> python -m environment_models.boot lattice --id lattice`

      This creates the Environment agent, waiting for Cell simulations to register.

      **NOTE:** If you didn't open the browser-based visualization, you can have the
      Environment agent open a "microscope" view onto the plate by launching it like this:

      `> ENVIRONMENT_ANIMATION=1 python -m environment_models.boot lattice --id lattice`

   2. Now start a Cell agent in a new tab:

      `> python -m environment_models.boot ecoli --id 1 --outer-id lattice`

      **Optional:** Supply additional arguments to set a variant, seed, and so on.
      Use the `-h` argument for help. 

5. Start as many cells as desired, each with its own unique id (agent name), and each in a
separate terminal tab.
You will see a message sent from the newly initialized simulation on the `environment_listen` topic:

   `<-- environment_listen: {'event': 'SIMULATION_INITIALIZED', 'id': '1'}`

   and in the environment tab you will see it has received a message from the new Cell simulation:

   `--> environment_listen: {u'event': u'SIMULATION_INITIALIZED', u'id': u'1'}`

6. Finally, run this in a separate "command" tab to start the simulation clock:

   `> python -m environment_models.boot trigger --id lattice`

7. To stop the simulation, run `shutdown` in the command tab:

   `> python -m environment_models.boot shutdown --id lattice`

## Agent Shepherd

An alternate way to start the simulation is to use the agent shepherd, which will manage the spawning and removal of agents with multiprocessing rather than launching each in its own tab. To do this first start the agent shepherd:

   `> python -m environment_models.boot shepherd`

Now that it is running you can start an experiment:

   `> python -m environment_models.boot experiment --number 3`

This will send four `ADD_AGENT` messages to the shepherd, one for the environment agent and three for the simulation agents. Note the `agent_id` for the lattice as you will need this for future control messages (like trigger and shutdown). These messages are received by the shepherd and you will see them all boot in the shepherd's tab. You still need to trigger execution, which requires the `agent_id` of the environment:

   `> python -m environment_models.boot trigger --id xxxxxx-xxxx-xxxxxxxxxx`

If you know that this is the only environment you are running, or you want to control all environments at once, you can omit the `--id` option:

   `> python -m environment_models.boot trigger`

Now that they are running, you can add new agents with `add`:

   `> python -m environment_models.boot add --id xxxxxx-xxxx-xxxxxxxxxx`

Or remove them with `remove` given an id. This can be just the prefix of the agent's id so you don't have to type the whole uuid:

   `> python -m environment_models.boot remove --id dgaf`

Finally, to shut down the experiment call `shutdown` as before:

   `> python -m environment_models.boot shutdown --id xxxxxx-xxxx-xxxxxxxxxx`

Notice this just shuts down the experiment, the shepherd is still running and a new experiment can be started. To shut down the shepherd process, just `Ctrl-C`.

## commands

To summarize, the list of agent control commands is:

* trigger - start the simulation
* pause - pause the simulation
* shutdown - shutdown the simulation
* divide - trigger cell division in the given simulation

Spawning specific agents:

* inner - start an inner agent
* outer - start an outer agent
* ecoli - start an ecoli agent
* lattice - start a two dimensional lattice agent

Shepherd oriented commands:

* shepherd - start the agent shepherd
* add - add an agent to the agent shepherd
* remove - remove an agent from the agent shepherd
* experiment - spawn an experiment from the agent shepherd
