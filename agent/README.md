# Agent

Distributed simulation of whole cell agents relative to a shared environment.

## Setup

The simulation is written in [Python 2.7.15](https://www.python.org/), and depends on [Kafka](https://kafka.apache.org/) for mediating communication between the different processes.

Kafka is a message passing system that allows decoupling of message senders and message receivers. It does this by providing two abstractions, a Consumer and a Producer. A Consumer can subscribe to any number of "topics" it will receive messages on, and a Producer can send to any topic it wishes. The topics provide a means of communication between processes that otherwise do not need to know any details of who is sending and receiving these messages.

Kafka is built on Zookeeper, which is a platform for coordinating access to an arbitrarily nested hierarchy of paths (called "nodes"). These are both meant to be deployed as a cluster on persistent servers somewhere, but there is a binary distribution that can be run locally. You will either have to have both of these running somewhere or access to a remote cluster running them in order to run this distributed simulation.

If you have access to a remote Kafka cluster, you can just specify the host as an option to the various boot scripts supplied here:

   `python -m agent.boot --host ip.to.remote.cluster:9092`

If you don't have access to a remote cluster, you can install Kafka locally.
 
1. If you don't already have it, [download Open JDK 8](https://jdk.java.net/8/).
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

   **Tip:** Use iTerm split windows to keep the Zookeeper and Kafka shells together.

   **Tip:** Create a shell alias like this in your bash profile:

   `alias zookeeper="cd ~/dev/kafka_2.11-2.0.0 && ./bin/zookeeper-server-start.sh ./config/zookeeper.properties"`

4. In another shell tab, in the same untar directory, start the Kafka server:

    `./bin/kafka-server-start.sh ./config/server.properties`

   **Tip:** Create a shell alias like this in your bash profile:

   `alias kafka="cd ~/dev/kafka_2.11-2.0.0 && ./bin/kafka-server-start.sh ./config/server.properties"`

You _can_ run stub agents using the command line tools `agent.boot` and `agent.control`, but
in practice we use the commands in [environment/README.md](../environment/README.md) to run
E. coli cell agents and their environment.
