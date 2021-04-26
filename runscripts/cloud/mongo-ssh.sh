#!/usr/bin/env bash
# ssh to the MongoDB server machine in GCE and tunnel to its MongoDB server port.
#
# COMMAND LINE ARGUMENTS:
#   ARG1 (optional): "bg" to run the ssh command in a background process.
#   ARG2 (optional): The MongoDB host name in GCE to ssh to.
#
# For port forwarding, this uses an explicit IPv4 local address 127.0.0.1 so if
# the port is in use ssh will fail and exit rather than just warning about it
# and using an IPv6 local address instead.

set -eu

HOST=${2:-mongo2}
ZONE=us-west1-b
PORT=27017
TUNNEL=127.0.0.1:$PORT:localhost:$PORT

if [ "${1-}" == bg ]
then
    gcloud compute ssh "$HOST" --zone="$ZONE" -- \
        -o ExitOnForwardFailure=yes -L $TUNNEL -nNT &
    ssh_pid="$!"

    echo "ssh port forwarding to $HOST in the background process: pid ${ssh_pid}"
    echo "Do 'kill ${ssh_pid}' to stop it."
else
    gcloud compute ssh "$HOST" --zone="$ZONE" -- \
        -o ExitOnForwardFailure=yes -L $TUNNEL
fi
