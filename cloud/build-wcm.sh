#!/bin/sh
# Use Google Cloud Build servers to build a personalized wcm-code Docker
# image and store it in the Google Container Registry.
#
# COMMAND LINE ARGUMENTS:
#   ARG1: Distinguishing ID prefix for the "wcm-code" tag; defaults to "$USER";
#     can identify a PR build, nightly build, git hash, ...
#   ARG2: Docker tag for the wcm-runtime image in GCR to build FROM; defaults to
#     "${ID}-wcm-runtime".
#     The named Docker image must already exist in the GCR project.
#
# ASSUMES: The current working dir is the wcEcoli/ project root.

set -eu

ID="${1:-$USER}"

WCM_RUNTIME="${2:-${ID}-wcm-runtime}"
WCM_CODE="${ID}-wcm-code"
GIT_HASH=$(git rev-parse HEAD)
GIT_BRANCH=$(git symbolic-ref --short HEAD)
TIMESTAMP=$(date '+%Y%m%d.%H%M%S')

echo "=== Cloud-building WCM code Docker Image ${WCM_CODE} on ${WCM_RUNTIME} ==="
echo "=== git hash ${GIT_HASH}, git branch ${GIT_BRANCH} ==="

# This needs a config file to identify the project files to upload and the
# Dockerfile to run.
gcloud builds submit --timeout=15m --config config-build-2-wcm-code.json \
    --substitutions="_WCM_RUNTIME=${WCM_RUNTIME},_WCM_CODE=${WCM_CODE},_GIT_HASH=${GIT_HASH},_GIT_BRANCH=${GIT_BRANCH},_TIMESTAMP=${TIMESTAMP}"
