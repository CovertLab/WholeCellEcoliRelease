#!/usr/bin/env bash
# Launch Sisyphus worker nodes with the given names on Google Compute Engine,
# provisioned with the capacity and service access needed for wcEcoli.

set -eu

DIR=$(dirname "$0")
python "${DIR}/gce_vms.py" --sisyphus --workflow "${WORKFLOW:?}" "$@"
