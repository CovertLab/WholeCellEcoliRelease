#!/bin/sh
# Use Google Cloud Build servers to build a personalized wcm-code Docker
# image and store it in the Google Container Registry.
#
# INPUTS:
#   ARG1: Distinguishing ID prefix for the "wcm-code" tag; defaults to "$USER";
#     can identify a PR build, nightly build, git hash, ...
#   ARG2: Docker tag for the wcm-runtime image in GCR to build FROM; defaults to
#     "wcm-runtime". This lets you test a new runtime that has new packages.
#     The named Docker image must already exist in the GCR project.
#
# ASSUMES: The current working dir is the wcEcoli/ project root.

set -eu

ID="${1:-$USER}"

WCM_RUNTIME="${2:-wcm-runtime}"
WCM_CODE="${ID}-wcm-code"

echo "=== Cloud-building WCM code ${WCM_CODE} on ${WCM_RUNTIME} ==="

# This needs a config file to identify the project files to upload and the
# Dockerfile to run.
gcloud builds submit --timeout=15m --config config-build-2-wcm-code.json \
    --substitutions="_WCM_RUNTIME=${WCM_RUNTIME},_WCM_CODE=${WCM_CODE}"
