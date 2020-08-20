#! /usr/bin/env bash
# Build the WCM Docker container images locally.
#
# ASSUMES: The current working dir is the wcEcoli/ project root.

set -eu

# On macOS, build with `NO_AVX2=1` to avoid the OpenBLAS 0.3.6+ self-test failures
# and bad results when building with Docker Desktop on macOS.
if [ "$(uname -s)" == Darwin ]; then NO_AVX2=1; else NO_AVX2=0; fi

# Docker image #1: The Python runtime environment.
docker build -f cloud/docker/runtime/Dockerfile -t wcm-runtime \
  --build-arg NO_AVX2=$NO_AVX2 .

# Docker image #2: The Whole Cell Model code on the runtime environment.
# See this Dockerfile for usage instructions.
docker build -f cloud/docker/wholecell/Dockerfile -t wcm-code .

# Docker image #3: The Whole Cell Model code with parameters calculated, ready for sims.
# See this Dockerfile for usage instructions.
docker build -f cloud/docker/full/Dockerfile -t wcm-full .
