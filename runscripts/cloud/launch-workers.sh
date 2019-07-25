#!/usr/bin/env bash
# Launch Sisyphus worker nodes with the given names on Google Compute Engine,
# provisioned with the capacity and service access needed for wcEcoli.

set -eu

NAMES=$@
PROJECT=allen-discovery-center-mcovert

# Could add --async option to return quickly if we needn't wait on successful creation.

gcloud compute \
       --project=$PROJECT \
       instances create $NAMES \
       --zone=us-west1-b \
       --machine-type=n1-standard-2 \
       --subnet=default \
       --network-tier=PREMIUM \
       --maintenance-policy=MIGRATE \
       --service-account=441871726775-compute@developer.gserviceaccount.com \
       --scopes=https://www.googleapis.com/auth/devstorage.read_only,https://www.googleapis.com/auth/logging.write,https://www.googleapis.com/auth/monitoring.write,https://www.googleapis.com/auth/servicecontrol,https://www.googleapis.com/auth/service.management.readonly,https://www.googleapis.com/auth/trace.append \
       --image-family=sisyphus-worker \
       --image-project=$PROJECT \
       --boot-disk-size=200GB \
       --boot-disk-type=pd-standard \
       --labels=role=sisyphus \
       --description='sisyphus worker'
