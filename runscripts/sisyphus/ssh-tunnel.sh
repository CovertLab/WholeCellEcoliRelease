#!/usr/bin/env bash
# ssh to the Gaia server and open tunnels to the Gaia and Kafka servers.
#
# REQUIRES: Your /etc/hosts needs the following line (not as a comment) so when
# the Kafka server reports its broker address as 'zookeeper-prime', connections
# to zookeeper-prime:9092 will also get forwarded to that host:
#     127.0.0.1   zookeeper-prime
#
# TODO(jerry): Configure ssh to handle zookeeper-prime w/o /etc/hosts?
# TODO(jerry): Add `-f`? `-o ExitOnForwardFailure=yes`?

set -eu

gcloud compute ssh gaia-base --zone=us-west1-b -- -L 24442:localhost:24442 -L 9092:zookeeper-prime:9092
