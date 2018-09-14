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

      `> python -m environment.boot lattice`

      This creates the Environment agent, waiting for Cell simulations to register.

      **NOTE:** If you didn't open the browser-based visualization, you can have the
      Environment agent open a "microscope" view onto the plate by launching it like this:

      `> ENVIRONMENT_ANIMATION=1 python -m environment.boot lattice`

   2. Now start a Cell agent in a new tab:

      `> python -m environment.boot ecoli --id 1`

      **Optional:** Supply additional arguments to set a variant, seed, and so on.
      Use the `-h` argument for help. 

5. Start as many cells as desired, each with its own unique id (agent name), and each in a
separate terminal tab.
You will see a message sent from the newly initialized simulation on the `environment_listen` topic:

   `<-- environment_listen: {'event': 'SIMULATION_INITIALIZED', 'id': '1'}`

   and in the environment tab you will see it has received a message from the new Cell simulation:

   `--> environment_listen: {u'event': u'SIMULATION_INITIALIZED', u'id': u'1'}`

6. Finally, run this in a separate "command" tab to start the simulation clock:

   `> python -m agent.boot trigger`

7. To stop the simulation, run `shutdown` in the command tab:

   `> python -m agent.boot shutdown`
