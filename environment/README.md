# Environment

Models of environments, which can be interfaced with cells via the agent module.
Includes command line tools `environment.boot` for running environment and cell agents
and `environment.control` for sending commands to agents.

## Setup

See the top-level [README.md](../README.md) for general setup instructions, and the
[agent README.md](../agent/README.md) for multi-agent simulation setup.

## Usage

1. To run Whole Cell E.coli simulations, you need to have the sim_data files. You can generate them via the
runParca manual runscript. In the wcEcoli directory:

    `> PYTHONPATH="$PWD" python runscripts/manual/runParca.py`

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

4. You can run an Environment agent and Cell agents directly from the command line although we usually do it via an agent Shepherd (see below).

   (**Tip:** Run each process in a new terminal tab. Use iTerm split windows to make it easy to watch them all at once.)

   In the first tab start an Environment model agent:

      `> python -m environment.boot --type lattice --id lattice`

      This creates the Environment agent, waiting for Cell simulations to register.
      You can optionally pass in a JSON `--config '{...}'` dictionary.

      **VARIATION:** If you didn't open the browser-based visualization, you can have the
      Environment agent open a "microscope" view onto the plate by launching it like this:

      `> ENVIRONMENT_ANIMATION=1 python -m environment.boot --type lattice --id lattice`

5. Now start a Cell agent in a new tab:

    `> python -m environment.boot --type ecoli --id 1 --outer-id lattice`

   You'll see messages like this one from the Cell agent to the Environment agent,
   declaring itself to the environment and giving its current state:

   `<-- environment-receive (CELL_DECLARE) [1]: {"inner_id": "1", "agent_config": {...}, "state": {"volume": 1.2, "environment_change": {}}, "event": "CELL_DECLARE", "agent_id": "lattice"}`

   In turn, the cell will receive messages like this:

   `--> cell-receive (ENVIRONMENT_SYNCHRONIZE) [1]: {u'inner_id': u'1', u'state': {u'time': 0}, u'outer_id': u'lattice', u'event': u'ENVIRONMENT_SYNCHRONIZE'}`

6. Start as many cells as desired, each with its own unique id (agent name), and each in a
separate terminal tab.

7. Finally, run `trigger` in a separate "command" tab to start the simulation clock:

   `> python -m environment.control trigger --id lattice`

8. To stop the simulation, run `shutdown` in the command tab:

   `> python -m environment.control shutdown --id lattice`

## Agent Shepherd

The current way to start the simulation is to use the agent Shepherd, which will manage the spawning and removal of agents as subprocesses rather than launching each in its own tab.

Clone the [CovertLab/shepherd](https://github.com/CovertLab/shepherd) repo and run:

   `> lein run`

Now that it is running you can start an experiment in another terminal tab:

   `> python -m environment.control experiment --number 3`

This will send four `ADD_AGENT` messages to the shepherd, one for the environment agent and three for the simulation agents. Note the `agent_id` for the lattice as you will need this for future control messages (like trigger and shutdown). These messages are received by the shepherd and you will see them all boot in the shepherd's tab. You still need to trigger execution, which requires the `agent_id` of the environment:

   `> python -m environment.control trigger --id xxxxxx-xxxx-xxxxxxxxxx`

If you know that this is the only environment you are running, or you want to control all environments at once, you can omit the `--id` option:

   `> python -m environment.control trigger`

Now that they are running, you can add new agents with `add`:

   `> python -m environment.control add --id xxxxxx-xxxx-xxxxxxxxxx`

Or remove them with `remove` given an id. This can be just the prefix of the agent's id so you don't have to type the whole uuid:

   `> python -m environment.control remove --id dgaf`

Finally, to shut down the experiment call `shutdown` as before:

   `> python -m environment.control shutdown --id xxxxxx-xxxx-xxxxxxxxxx`

Notice this just shuts down the experiment, the shepherd is still running and a new experiment can be started. To shut down the shepherd process, just `Ctrl-C`.

## command summary

This is just a summary.
Use the `-h` argument to get complete usage help on these command line programs.

The environment.boot commands run an agent in the current shell tab:

* ecoli - an ecoli cell agent
* lattice - a two dimensional lattice environment agent
* chemotaxis - a chemotaxis surrogate that can move up glucose gradients within a chemotaxis_experiment

The environment.control commands include:

* trigger - start/resume the simulation clock
* pause - pause the simulation clock
* shutdown - shutdown the simulation

Some environment.control commands require an [agent shepherd](https://github.com/CovertLab/shepherd), including:

* add - ask the shepherd to spawn an agent and add it to an environment
* remove - ask the shepherd to remove an agent
* experiment - ask the shepherd to spawn a lattice and multiple cell agents
