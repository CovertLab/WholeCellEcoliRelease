#!/usr/bin/env bash
# ssh to the Gaia server and open a TCP tunnel to the Gaia server.
#
# TODO(jerry): Add `-f`? `-o ExitOnForwardFailure=yes`?

set -eu

gcloud compute ssh gaia-base --zone=us-west1-b -- -L 24442:localhost:24442
