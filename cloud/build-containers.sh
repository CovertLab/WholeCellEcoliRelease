#!/bin/sh
# Build the WCM Docker container images locally.
#
# ASSUMES: The current working dir is the wcEcoli/ project root.

set -eu

# Docker image #1: The Python runtime environment.
docker build -f cloud/docker/runtime/Dockerfile -t wcm-runtime .

# Docker image #2: The Whole Cell Model code on the runtime environment.
# See this Dockerfile for usage instructions.
docker build -f cloud/docker/wholecell/Dockerfile -t wcm-code .
