#!/usr/bin/env bash
# Build the WCM Docker container images locally.
#
# ASSUMES: The current working dir is the wcEcoli/ project root.
#
# Add the `docker build` option `--build-arg from=ABC` to name a different "FROM" image.

set -eu

# To avoid OpenBLAS self-test failures when compiling it in Docker Desktop on
# macOS, build the Image with `COMPILE_BLAS=0` (=> don't compile OpenBLAS) or
# with `NO_AVX2=1` (=> compile it to not use AVX2 vector instructions).
COMPILE_BLAS=0
if [ "$(uname -s)" == Darwin ]; then NO_AVX2=1; else NO_AVX2=0; fi
echo COMPILE_BLAS=$COMPILE_BLAS, NO_AVX2=$NO_AVX2

WCM_RUNTIME=${USER}-wcm-runtime
WCM_CODE=${USER}-wcm-code
GIT_HASH=$(git rev-parse HEAD)
GIT_BRANCH=$(git symbolic-ref --short HEAD)
TIMESTAMP=$(date '+%Y%m%d.%H%M%S')

# Docker image #1: The Python runtime environment.
docker build -f cloud/docker/runtime/Dockerfile -t "${WCM_RUNTIME}" \
  --build-arg COMPILE_BLAS=$COMPILE_BLAS \
  --build-arg NO_AVX2=$NO_AVX2 \
  .

# Docker image #2: The Whole Cell Model code on the runtime environment.
# See this Dockerfile for usage instructions.
docker build -f cloud/docker/wholecell/Dockerfile -t "${WCM_CODE}" \
  --build-arg from="${WCM_RUNTIME}" \
  --build-arg git_hash="${GIT_HASH}" \
  --build-arg git_branch="${GIT_BRANCH}" \
  --build-arg timestamp="${TIMESTAMP}" .
