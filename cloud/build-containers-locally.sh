#!/usr/bin/env bash
# Build the WCM Docker container images locally.
#
# ASSUMES: The current working dir is the wcEcoli/ project root.
#
# Add the `docker build` option `--build-arg from=ABC` to name a different "FROM" image.

set -eu

# Build with `NO_AVX2=1` on macOS to avoid the OpenBLAS 0.3.6+ fatal self-test
# errors and runtime bad results in Docker Desktop on macOS. This workaround
# won't help if you build an image on Linux then run it on macOS. See the Dockerfile.
if [ "$(uname -s)" == Darwin ]; then NO_AVX2=1; else NO_AVX2=0; fi
echo NO_AVX2=$NO_AVX2

GIT_HASH=$(git rev-parse HEAD)
GIT_BRANCH=$(git symbolic-ref --short HEAD)
TIMESTAMP=$(date '+%Y%m%d.%H%M%S')

# Docker image #1: The Python runtime environment.
docker build -f cloud/docker/runtime/Dockerfile -t wcm-runtime \
  --build-arg NO_AVX2=$NO_AVX2 .

# Docker image #2: The Whole Cell Model code on the runtime environment.
# See this Dockerfile for usage instructions.
docker build -f cloud/docker/wholecell/Dockerfile -t wcm-code \
  --build-arg git_hash="${GIT_HASH}" \
  --build-arg git_branch="${GIT_BRANCH}" \
  --build-arg timestamp="${TIMESTAMP}" .
