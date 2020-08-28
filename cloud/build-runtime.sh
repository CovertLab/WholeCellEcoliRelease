#!/bin/sh
# Use Google Cloud Build servers to build a personalized wcm-runtime Docker
# image and store it in the Google Container Registry.
#
# COMMAND LINE ARGUMENTS:
#   ARG1: Docker tag for the wcm-runtime image in GCR to build; defaults to
#     "${USER}-wcm-runtime".
#
# ASSUMES: The current working dir is the wcEcoli/ project root.

set -eu

PROJECT="$(gcloud config get-value core/project)"
WCM_RUNTIME="${1:-${USER}-wcm-runtime}"
TAG="gcr.io/${PROJECT}/${WCM_RUNTIME}"

echo "=== Cloud-building WCM runtime Docker Image: ${TAG} ==="

# This needs only one payload file so copy it in rather than using a config at
# the project root which would upload the entire project.
cp requirements.txt cloud/docker/runtime/
gcloud builds submit --timeout=2h --tag "${TAG}" cloud/docker/runtime/
rm cloud/docker/runtime/*requirements.txt
