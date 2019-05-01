# Agent

Distributed simulation of whole cell agents relative to a shared environment.

## Setup

The simulation is written in [Python 2.7](https://www.python.org/), and depends on [Kafka](https://kafka.apache.org/) for mediating communication between the different processes.

Kafka is a message passing system that allows decoupling of message senders and message receivers. It does this by providing two abstractions, a Consumer and a Producer. A Consumer can subscribe to any number of "topics" it will receive messages on, and a Producer can send to any topics it wishes. Topics are communication "channels" between processes that otherwise do not need to know who is sending and receiving these messages.

Kafka is built on Zookeeper, which is a service for synchronizing access to a hierarchy of key/value pairs called "nodes". The Kafka and Zookeeper services are meant to be deployed in multiple instances in a server cluster, but there is a binary distribution that you can run locally for development.

If you have access to a remote Kafka service, just specify the host as an option to the various boot scripts supplied here:

   `python -m agent.boot --host ip.to.remote.cluster:9092`

If you don't have access to a remote Kafka service, run Kafka and Zookeeper servers locally.

1. If you don't already have it, [download Open JDK 8](https://jdk.java.net/8/).
   1. Run the JDK installer.
   2. Set `JAVA_HOME` in your shell setup file (e.g. `.bash_profile` or `.profile`):

      `export JAVA_HOME=$(/usr/libexec/java_home)`

   3. Restart your shell to get the `JAVA_HOME` setting.
   4. Test it by running:

      `java -version`

      That should print something like

      `java version "1.8.0_202-ea"`

2. [Download the Apache Kafka server software](https://www.apache.org/dyn/closer.cgi?path=/kafka/2.0.0/kafka_2.11-2.0.0.tgz).
   1. **Optional:** Test the integrity of the downloaded `.tgz` file by [computing its SHA-1 checksum](https://www.apache.org/info/verification.html)
and pasting the checksum into [the verification page](https://www.apache.org/info/verification.html).
   2. Untar the `.tgz` file to a suitable directory such as `~/dev/kafka_2.11-2.0.0/`:

      `tar xvf kafka_2.11-2.0.0.tgz`

3. In your untar directory, start Zookeeper first:

   `bin/zookeeper-server-start.sh config/zookeeper.properties`

   Zookeeper will keep running until forced to stop.

   **Tip:** Remember to shut down Kafka before Zookeeper.

   **Tip:** Use iTerm split windows to keep the Zookeeper and Kafka shells together.

   **Tip:** Create a shell alias like this in your bash profile:

   `alias zookeeper="cd ~/dev/kafka_2.11-2.0.0 && bin/zookeeper-server-start.sh config/zookeeper.properties"`

4. In another shell tab, in the same untar directory, start the Kafka server:

   `bin/kafka-server-start.sh config/server.properties --override listeners=PLAINTEXT://127.0.0.1:9092`

   **Note:** Overriding the "listeners" address like this allows connections to the Kafka server to withstand network DHCP address changes and the like.

   **Tip:** Create a shell alias like this in your bash profile:

   `alias kafka="cd ~/dev/kafka_2.11-2.0.0 && bin/kafka-server-start.sh config/server.properties --override listeners=PLAINTEXT://127.0.0.1:9092"`

   Kafka will keep running until forced to stop.
   If it has no clients, `Control-C` usually works to stop it.
   If it refuses to stop, use a `kill -9` command to do it.

You _can_ run stub agents using the command line tools `agent.boot` and `agent.control`, but
in practice we use the commands in [environment/README.md](../environment/README.md) to run
_E. coli_ cell agents and their environment.
