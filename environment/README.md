# Environment

Distributed simulation of whole cell agents relative to a shared environment.

## Setup

The simulation is written in [Python 2.7.15](https://www.python.org/), and depends on [Kafka](https://kafka.apache.org/) for mediating communication between the different processes.

To install Kafka you can download the packages [here](https://www.apache.org/dyn/closer.cgi?path=/kafka/2.0.0/kafka_2.11-2.0.0.tgz). 

Once untarred, start Zookeeper first:

    # in your untar directory
    > ./bin/zookeeper-server-start.sh ./config/zookeeper.properties

then start Kafka:

    > ./bin/kafka-server-start.sh ./config/server.properties

With this you should be ready to go.

## Usage

You will want some way to run multiple command line processes, as each component of the system claims the execution thread. This can be in multiple terminal tabs locally, or running on multiple VM's in a cloud environment. Here I will refer to them as tabs.

The command line interface for the system is is `environment/boot.py`, and has a number of options depending on the role you are interacting with.

The available roles are

* `outer` - the larger environmental context
* `inner` - each individual simulation
* `trigger` - start the execution of the system
* `shutdown` - stop the execution of the system

In the first tab start the environmental context:

    0> python environment/boot.py outer

This has started the environmental process and is waiting for simulations initialize and register. Let's do that now. In a new tab:

    1> python environment/boot.py inner --id 1

When starting individual simulations the `id` argument is required. Assigning a unique id to each simulation allows the environment to communicate in specific ways with each simulation. You will see a message sent from the newly initialized simulation on the `environment_listen` topic:

    <-- environment_listen: {'event': 'SIMULATION_INITIALIZED', 'id': '1'}

and also if you go back to the environment tab you will see it has received a message from the new simulation:

    --> environment_listen: {u'event': u'SIMULATION_INITIALIZED', u'id': u'1'}

Let's start another one in another tab:

    2> python environment/boot.py inner --id 2
    <-- environment_listen: {'event': 'SIMULATION_INITIALIZED', 'id': '2'}

We will see this message has also reached the environmental process:

    --> environment_listen: {u'event': u'SIMULATION_INITIALIZED', u'id': u'2'}

Now that we have a couple of simulations running, let's start the execution. In yet another tab (the `control` tab):

    3> python environment/boot.py trigger

You will see this sends a message to the `environment_control` topic:

    <-- environment_control: {'event': 'TRIGGER_EXECUTION'}

which is then received by the environment:

    --> environment_control: {u'event': u'TRIGGER_EXECUTION'}

Now things have started to run. You will see in the environment tab that it has sent two messages out along the `environment_broadcast` topic:

    <-- environment_broadcast: {'run_for': 1, 'molecule_ids': ['blue', 'green', 'red', 'yellow'], 'id': u'1', 'message_id': 0, 'concentrations': {'blue': 12, 'green': 11, 'red': 44, 'yellow': 5}, 'event': 'ENVIRONMENT_UPDATED'}
    <-- environment_broadcast: {'run_for': 1, 'molecule_ids': ['blue', 'green', 'red', 'yellow'], 'id': u'2', 'message_id': 0, 'concentrations': {'blue': 12, 'green': 11, 'red': 44, 'yellow': 5}, 'event': 'ENVIRONMENT_UPDATED'}

These are in turn received by the simulation, each responding to the message given by its id:

    --> environment_broadcast: {u'run_for': 1, u'molecule_ids': [u'blue', u'green', u'red', u'yellow'], u'id': u'1', u'message_id': 0, u'concentrations': {u'blue': 12, u'green': 11, u'yellow': 5, u'red': 44}, u'event': u'ENVIRONMENT_UPDATED'}

This message tells the simulation what its molecule ids are, the respective concentrations, and how long to run for. Upon receiving this message it begins the simulation and runs it for the given amount of time (`run_for`). Once it completes this segment of execution, it responds with a message on the `environment_listen` topic detailing the local changes to concentrations it has calculated:

    <-- environment_listen: {'message_id': 0, 'time': 1, 'changes': {u'blue': 1, u'green': 6, u'yellow': 3, u'red': 2}, 'event': 'SIMULATION_ENVIRONMENT', 'id': '1'}

The environment will wait until it has received messages from all its registered simulations, at which point it will perform the work to integrate the separate changes in concentrations it has received and respond with a new message to each simulation with its now updated environment. These are of the same form as the 'ENVIRONMENT_UPDATED' messages above.

Once this has gone on awhile and you are ready to stop the simulation, you can call shutdown from the command tab:

    3> python environment/boot.py shutdown

This sends an `ENVIRONMENT_SHUTDOWN` message on the `environment_control` topic:

    <-- environment_control: {'event': 'SHUTDOWN_ENVIRONMENT'}

The environment receives this message and waits until it receives all outstanding messages from the individual simulations, and then sends a shutdown message to each one:

    --> environment_control: {u'event': u'SHUTDOWN_ENVIRONMENT'}
    <-- environment_broadcast: {'id': u'1', 'event': 'SHUTDOWN_SIMULATION'}
    <-- environment_broadcast: {'id': u'2', 'event': 'SHUTDOWN_SIMULATION'}

At this point all three processes have exited.
