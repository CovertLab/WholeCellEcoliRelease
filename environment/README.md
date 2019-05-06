# Environment

Models of environments, which can be interfaced with cells via the agent module.
Includes command line tools `environment.boot` for running environment and cell agents
and `environment.control` for sending commands to agents.

## Setup

See the top-level [README.md](../README.md) for general setup instructions, and the
[agent README.md](../agent/README.md) for multi-agent simulation setup.

1. To run Whole Cell E.coli simulations, you need to have the sim_data files. You can generate them via the
runFitter manual runscript. In the wcEcoli directory:

    `> PYTHONPATH="$PWD" python runscripts/manual/runFitter.py`

## Zookeeper and Kafka

2. See [agent/README.md](../agent/README.md) for instructions to set up, start, and stop your Zookeeper and Kafka servers. To recap:

   1. Start Zookeeper in the directory where you untarred the Kafka and Zookeeper software:

      `> bin/zookeeper-server-start.sh config/zookeeper.properties`

   2. Then start the Kafka server in another shell tab in the same directory:

      `> bin/kafka-server-start.sh config/server.properties --override listeners=PLAINTEXT://127.0.0.1:9092`

3. **Optional:** Start the "Lens" environment visualization server and browser window per the instructions on the [CovertLab/shepherd](https://github.com/CovertLab/shepherd) page. To recap:

   1. Run the Lens visualization server in the root directory of the [CovertLab/shepherd](https://github.com/CovertLab/shepherd) repository:

      `> lein run -m shepherd.lens`

   2. Open a browser window onto [http://localhost:33332](http://localhost:33332)

## One Agent Per Terminal Tab (esp. for debugging)

**Note:** This section describes how to launch Environment and Cell agent processes directly
from the command line. We usually use the "Agent Shepherd" approach described in the next section,
but the direct approach lets you launch agents under a debugger. You can use it for one or more
agents together with others running under a shepherd.

**Tip:** Run each process in a new terminal tab. Use iTerm split windows to make it easy to watch them all at once.

4. In the first terminal tab, launch an Environment agent:

      `> python -m environment.boot --type lattice --id lattice`

      The Environment agent will wait for Cell simulation agents to register.
      You can optionally pass in a JSON `--config '{...}'` dictionary.

5. Now start a Cell agent in a new tab:

   `> python -m environment.boot --type ecoli --id 1 --outer-id lattice`

   Vary the agent type and other parameters as needed. Each agent needs an `id` that's unique among the
   currently running agents.

   You'll see messages like this one from the Cell agent to the Environment agent,
   declaring itself to the environment and giving its current state:

   `<-- environment-receive (CELL_DECLARE) [1]: {"inner_id": "1", "agent_config": {...}, "state": {"volume": 1.2, "environment_change": {}}, "event": "CELL_DECLARE", "agent_id": "lattice"}`

   In turn, the cell will receive messages like this:

   `--> cell-receive (ENVIRONMENT_SYNCHRONIZE) [1]: {u'inner_id': u'1', u'state': {u'time': 0}, u'outer_id': u'lattice', u'event': u'ENVIRONMENT_SYNCHRONIZE'}`

6. To debug an agent using the PyCharm debugger, first create a Run/Debug Configuration with the
module name and parameters on the command lines above, e.g.

   ```
   Name: cell agent
   
   Module name: environment.boot
   Parameters: ecoli --id 40
   Python interpreter: Python 2.7 (wcEcoli2)
   Working directory: /Users/YOU/dev/wcEcoli
   ☑︎ Add content roots to PYTHONPATH
   ☑︎ Add source roots to PYTHONPATH
   ☑︎ Run with Python console
   ```

   Then set your breakpoints and invoke the `Debug 'cell agent'` menu command.

6. Start as many cells as desired, each with its own unique `id` agent name; each in a
separate terminal tab or PyCharm run/debug tab.

7. Finally, use another terminal tab to start the simulation running:

   `> python -m environment.control run`

   You can `pause` and `run` it whenever you want.

8. To shut down the simulation, run `shutdown` in the command tab:

   `> python -m environment.control shutdown`

## Agent Shepherd

The usual way to start the simulation is to use the agent "Shepherd", which is a process
that spawns agents in subprocesses (as requested via Kafka messages) so you don't have to
launch each agent in its own terminal tab.
Furthermore, this enables cell division wherein a cell agent process ends and two
new ones begin.
But to debug an agent, see the "One Agent Per Terminal Tab" instructions, above.

Clone the [CovertLab/shepherd](https://github.com/CovertLab/shepherd) repo and run:

   `> lein run`

Use a "Lens" browser page to view the agents in action. To do this, open another shell
tab onto the shepherd repo directory and run:

   `> lein run -m shepherd.lens`

then open a browser window onto [http://localhost:33332/](http://localhost:33332/)

Now you can start a virtual microscope experiment in a "command" terminal tab:

   `> python -m environment.control experiment --number 3`

This will send four `ADD_AGENT` messages to the shepherd: one for the _lattice environment_ agent and three for the _cell simulation_ agents. Note the `agent_id` for the lattice as you will need this for future control messages (like `run` and `shutdown`). These messages are received by the shepherd and you will see all the agents' logs in the "shepherd" tab.

You can `run`/`pause` the simulation at will:

   `> python -m environment.control run`

   `> python -m environment.control pause`

You can add another cell agent:

   `> python -m environment.control add`

(If you're running multiple environment agents, you can specify a lattice environment agent id via the `--id` option.)

You can remove a cell agent using the prefix of the agent's id (you don't have to type the whole id):

   `> python -m environment.control remove --id dgaf`

Finally, to shut down the experiment, run `shutdown`:

   `> python -m environment.control shutdown`

Notice this just shuts down the agents. The Shepherd is still running and ready for a new experiment.
Use `Ctrl-C` to stop the Shepherd and Lens processes.

## command summary

This is just a summary.
Use the `-h` argument to get complete usage help on these command line programs.

The environment.boot commands run an agent in the current shell tab:

* ecoli - an ecoli cell agent
* lattice - a two dimensional lattice environment agent
* chemotaxis - a chemotaxis surrogate that can move up glucose gradients within a chemotaxis_experiment

The environment.control commands include:

* run - start/resume the simulation clock
* pause - pause the simulation clock
* shutdown - shutdown the simulation

Some environment.control commands require an [agent shepherd](https://github.com/CovertLab/shepherd), including:

* add - ask the shepherd to spawn an agent and add it to an environment
* remove - ask the shepherd to remove an agent
* experiment - ask the shepherd to spawn a lattice and multiple cell agents
