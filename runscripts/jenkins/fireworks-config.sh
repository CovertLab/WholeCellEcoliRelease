#! /usr/bin/env sh

# Create a new launchpad to allow fireworks to access a database to run a workflow.
# $WC_MONGO_USER, $WC_MONGO_PW, $WC_MONGO_CLUSTER should be set by Jenkins.
# Passing in a new database name as the first arg to this script will create
# a new database in the cluster for easy expansion for new Jenkins jobs.
# Be sure not to reuse a database name for another job or else jobs may wipe
# other workflows.

set -eu

WC_MONGO_DB=$1
LOG_DIR="/scratch/groups/mcovert/jenkins/fireworks/logs/launchpad"

mkdir -p $LOG_DIR
{
  echo "authsource: admin"
  echo "host: mongodb+srv://${WC_MONGO_USER}:${WC_MONGO_PW}@${WC_MONGO_CLUSTER}/${WC_MONGO_DB}?retryWrites=true&w=majority"
  echo "logdir: $LOG_DIR"
  echo "mongoclient_kwargs: {}"
  echo "name: null"
  echo "password: null"
  echo "port: null"
  echo "ssl: false"
  echo "ssl_ca_certs: null"
  echo "ssl_certfile: null"
  echo "ssl_keyfile: null"
  echo "ssl_pem_passphrase: null"
  echo "strm_lvl: INFO"
  echo "uri_mode: true"
  echo "user_indices: []"
  echo "username: null"
  echo "wf_user_indices: []"
} > my_launchpad.yaml
