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

### # Commented out until needed.
### # 3. The full WCM with Parca output, ready to run sims, esp. agent sims.
### # This build doesn't need to upload any payload files.
### gcloud builds submit --timeout=90m --tag gcr.io/allen-discovery-center-mcovert/wcm-full cloud/docker/full
