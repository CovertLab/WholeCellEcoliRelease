#!/bin/sh
# Use Google Cloud Build servers to build the tiered WCM Docker images and store
# them in the Google Container Registry.
#
# ASSUMES: The current working dir is the wcEcoli/ project root.

set -eu

# 1. The Python runtime environment with the default image tag.
cloud/build-runtime.sh

# 2. The Whole Cell Model code with the default user-specific tag, building
# FROM the default runtime.
cloud/build-wcm.sh
