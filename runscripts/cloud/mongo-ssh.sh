#!/usr/bin/env bash
# ssh to the MongoDB server machine and tunnel to the MongoDB server processes.
#
# For port forwarding, this uses an explicit local address 127.0.0.1:27017 so
# it will fail and exit if that port is in use rather than just warning and
# using an IPv6 address.
#
# TODO: Take an optional port argument in case the user is using local port
# 27017 for a local MongoDB server or something.

set -eu

gcloud compute ssh mongo-prime --zone=us-west1-b -- \
    -o ExitOnForwardFailure=yes -L 127.0.0.1:27017:localhost:27017
