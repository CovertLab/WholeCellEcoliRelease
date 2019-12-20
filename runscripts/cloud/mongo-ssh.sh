#!/usr/bin/env bash
# ssh to the MongoDB server machine and tunnel to the MongoDB server processes.
#
# TODO(jerry): Add `-f`? `-o ExitOnForwardFailure=yes`?

set -eu

gcloud compute ssh mongo-prime --zone=us-west1-b -- -L 27017:localhost:27017
