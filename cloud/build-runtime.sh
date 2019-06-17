#!/bin/sh
# Use Google Cloud Build servers to build a personalized wcm-runtime Docker
# image and store it in the Google Container Registry.
#
# INPUTS:
#   ARG1: Docker tag for the wcm-runtime image in GCR to build; defaults to
#     "wcm-runtime". This lets you build a test runtime that has new packages
#     without breaking anyone else.
#
# ASSUMES: The current working dir is the wcEcoli/ project root.

set -eu

WCM_RUNTIME="${1:-wcm-runtime}"

echo "=== Cloud-building WCM runtime: ${WCM_RUNTIME} ==="

# This needs one payload file so copy it in rather than using a config at the
# project root which would upload the entire project.
cp requirements.txt cloud/docker/runtime/
gcloud builds submit --timeout=2h \
    --tag "gcr.io/allen-discovery-center-mcovert/${WCM_RUNTIME}" \
    cloud/docker/runtime/
rm -f cloud/docker/runtime/requirements.txt
