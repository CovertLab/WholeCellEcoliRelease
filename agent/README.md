# Agent

Distributed simulation of whole cell agents relative to a shared environment.

## Setup

The simulation is written in [Python 2.7.15](https://www.python.org/), and depends on [Kafka](https://kafka.apache.org/) for mediating communication between the different processes.

Kafka is a message passing system that allows decoupling of message senders and message receivers. It does this by providing two abstractions, a Consumer and a Producer. A Consumer can subscribe to any number of "topics" it will receive messages on, and a Producer can send to any topic it wishes. The topics provide a means of communication between processes that otherwise do not need to know any details of who is sending and receiving these messages.

Kafka is built on Zookeeper, which is a platform for coordinating access to an arbitrarily nested hierarchy of paths (called "nodes"). These are both meant to be deployed as a cluster on persistent servers somewhere, but there is a binary distribution that can be run locally. You will either have to have both of these running somewhere or access to a remote cluster running them in order to run this distributed simulation.

If you have access to a remote Kafka cluster, you can just specify the host as an option to the various boot scripts supplied here:

    python -m agent.boot --host ip.to.remote.cluster:9092

If you don't have access to a remote cluster, you can install Kafka locally.
 
1. If you don't already have it, [download Java JDK 8](http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html).
   1. Run the JDK installer.
   2. Set `JAVA_HOME` in your shell setup file (e.g. `.bash_profile` or `.profile`):

      `export JAVA_HOME=$(/usr/libexec/java_home)`

   3. Restart your shell to get the `JAVA_HOME` setting.
   4. Test it

      `java -version`

      That should print something like

      `java version "1.8.0_181"`

2. [Download the Apache Kafka server software](https://www.apache.org/dyn/closer.cgi?path=/kafka/2.0.0/kafka_2.11-2.0.0.tgz).
   1. [Optional] Test the integrity of the downloaded `.tgz` file by [computing its SHA-1 checksum](https://www.apache.org/info/verification.html)
and pasting the checksum into [the verification page](https://www.apache.org/info/verification.html).
   2. Untar the `.tgz` file in a suitable directory.

      `tar xvf kafka_2.11-2.0.0.tgz`

3. In your untar directory, start Zookeeper first:

   `./bin/zookeeper-server-start.sh ./config/zookeeper.properties`

   It will keep running until forced to shut down.
4. In another shell tab, in the same untar directory, start the Kafka server:

    `./bin/kafka-server-start.sh ./config/server.properties`

With this you should be ready to go.

## Usage

You will want some way to run multiple command line processes, as each component of the system claims the execution thread. This can be in multiple terminal tabs locally, or running on multiple VM's in a cloud environment. Here I will refer to them as tabs.

The command line interface for the system is is `agent/boot.py`, and has a number of options depending on the role you are interacting with.

The available roles are

* `outer` - the larger environmental context
* `inner` - each individual simulation
* `trigger` - start the execution of the system
* `shutdown` - stop the execution of the system

In the first tab start the environmental context:

    0> python -m agent.boot outer

This has started the environmental process and is waiting for simulations initialize and register. Let's do that now. In a new tab:

    1> python -m agent.boot inner --id 1

When starting individual simulations the `id` argument is required. Assigning a unique id to each simulation allows the environment to communicate in specific ways with each simulation. You will see a message sent from the newly initialized simulation on the `environment_listen` topic:

    <-- environment_listen: {'event': 'SIMULATION_INITIALIZED', 'id': '1'}

and also if you go back to the environment tab you will see it has received a message from the new simulation:

    --> environment_listen: {u'event': u'SIMULATION_INITIALIZED', u'id': u'1'}

Let's start another one in another tab:

    2> python -m agent.boot inner --id 2
    <-- environment_listen: {'event': 'SIMULATION_INITIALIZED', 'id': '2'}

We will see this message has also reached the environmental process:

    --> environment_listen: {u'event': u'SIMULATION_INITIALIZED', u'id': u'2'}

Now that we have a couple of simulations running, let's start the execution. In yet another tab (the `control` tab):

    3> python -m agent.boot trigger

You will see this sends a message to the `environment_control` topic:

    <-- environment_control: {'event': 'TRIGGER_EXECUTION'}

which is then received by the environment:

    --> environment_control: {u'event': u'TRIGGER_EXECUTION'}

Now things have started to run. You will see in the environment tab that it has sent two messages out along the `environment_broadcast` topic:

    <-- environment_broadcast: {'run_for': 1, 'id': u'1', 'message_id': 0, 'concentrations': {'blue': 12, 'green': 11, 'red': 44, 'yellow': 5}, 'event': 'ENVIRONMENT_UPDATED'}
    <-- environment_broadcast: {'run_for': 1, 'id': u'2', 'message_id': 0, 'concentrations': {'blue': 12, 'green': 11, 'red': 44, 'yellow': 5}, 'event': 'ENVIRONMENT_UPDATED'}

These are in turn received by the simulation, each responding to the message given by its id:

    --> environment_broadcast: {u'run_for': 1, u'id': u'1', u'message_id': 0, u'concentrations': {u'blue': 12, u'green': 11, u'yellow': 5, u'red': 44}, u'event': u'ENVIRONMENT_UPDATED'}

This message tells the simulation what its molecule ids are, the respective concentrations, and how long to run until. Upon receiving this message it begins the simulation and runs until it reaches the given time step (`run_until`). Once it completes this segment of execution, it responds with a message on the `environment_listen` topic detailing the local changes to concentrations it has calculated:

    <-- environment_listen: {'message_id': 0, 'time': 1, 'changes': {u'blue': 1, u'green': 6, u'yellow': 3, u'red': 2}, 'event': 'SIMULATION_ENVIRONMENT', 'id': '1'}

The environment will wait until it has received messages from all its registered simulations, at which point it will perform the work to integrate the separate changes in concentrations it has received and respond with a new message to each simulation with its now updated environment. These are of the same form as the 'ENVIRONMENT_UPDATED' messages above.

Once this has gone on awhile and you are ready to stop the simulation, you can call shutdown from the command tab:

    3> python -m agent.boot shutdown

This sends an `ENVIRONMENT_SHUTDOWN` message on the `environment_control` topic:

    <-- environment_control: {'event': 'SHUTDOWN_ENVIRONMENT'}

The environment receives this message and waits until it receives all outstanding messages from the individual simulations, and then sends a shutdown message to each one:

    --> environment_control: {u'event': u'SHUTDOWN_ENVIRONMENT'}
    <-- environment_broadcast: {'id': u'1', 'event': 'SHUTDOWN_SIMULATION'}
    <-- environment_broadcast: {'id': u'2', 'event': 'SHUTDOWN_SIMULATION'}

At this point all three processes have exited.


**Tip:** You can shut down an individual agent (say #1) like this:

    python -m agent.boot shutdown --id 1

That's handy if the Outer agent raised an exception and exited.
