#!/usr/bin/env bash
# Build the WCM Docker container images locally.
#
# ASSUMES: The current working dir is the wcEcoli/ project root.

set -eu

# Docker image #1: The Python runtime environment.
cloud/locally-build-runtime.sh

# Docker image #2: The Whole Cell Model code on the runtime environment.
# See this Dockerfile for usage instructions.
cloud/locally-build-wcm.sh
