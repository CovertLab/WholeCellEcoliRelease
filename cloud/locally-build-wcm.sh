#!/usr/bin/env bash
# Build the wcm-code Docker container Image locally.
#
# ASSUMES: The current working dir is the wcEcoli/ project root.

set -eu

WCM_RUNTIME=${USER}-wcm-runtime
WCM_CODE=${USER}-wcm-code
GIT_HASH=$(git rev-parse HEAD)
GIT_BRANCH=$(git symbolic-ref --short HEAD)
TIMESTAMP=$(date '+%Y%m%d.%H%M%S')

mkdir -p source-info
git diff HEAD > source-info/git_diff.txt

# Docker image #2: The Whole Cell Model code on the runtime environment.
# See this Dockerfile for usage instructions.
docker build -f cloud/docker/wholecell/Dockerfile -t "${WCM_CODE}" \
  --build-arg from="${WCM_RUNTIME}" \
  --build-arg git_hash="${GIT_HASH}" \
  --build-arg git_branch="${GIT_BRANCH}" \
  --build-arg timestamp="${TIMESTAMP}" .

rm source-info/git_diff.txt
