#!/bin/sh
# Use Cloud Build servers to build the layered WCM Docker images and store them
# in the Container Registry.
#
# ASSUMES: The current working dir is the wcEcoli/ project root.
#
# TODO: Try using `gcloud builds submit --substitutions TAG_NAME=...` for
# different developers and CI builds to create images containing different code.
# That should work for the config files but we'll have to substitute into the
# Dockerfile "FROM" statements.

set -eu

# 1. The Python runtime environment.
# This needs one payload file so copy it in rather than using a config at the
# project root which would upload the entire project.
cp requirements.txt cloud/docker/runtime/
gcloud builds submit --timeout=2h --tag gcr.io/allen-discovery-center-mcovert/wcm-runtime cloud/docker/runtime/
rm cloud/docker/runtime/requirements.txt

# 2. The Whole Cell Model code.
# This needs a config file to identify the project files to upload and the
# Dockerfile to run.
gcloud builds submit --timeout=15m --config config-build-2-wcm-code.json

# 3. The full WCM with Parca output, ready to run sims.
# This build doesn't need to upload any payload files.
gcloud builds submit --timeout=90m --tag gcr.io/allen-discovery-center-mcovert/wcm-full cloud/docker/full
