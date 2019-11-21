# Container image #3: wcm-full.
# This Dockerfile builds a container image with the full Whole Cell Model, ready
# to run simulations, layered on the wcm-code image.
#
# To build locally from the wcEcoli/ project root directory:
#
#     > docker build -f cloud/docker/full/Dockerfile -t wcm-full .
#
# To build using the Cloud Build service and store in the Container Registry,
# from this directory:
#
#     > gcloud builds submit --timeout=90m --tag gcr.io/allen-discovery-center-mcovert/wcm-full .
#
# After building locally you can start up a new container from the image:
#
#     > docker run --name wcm --user "$(id -u):$(id -g)" -it --rm wcm-full
#
# It will start a shell where you can execute commands:
#
#     # pytest
#
# If this succeeds you should be good to go, e.g.:
#
#     # python runscripts/manual/runSim.py


FROM gcr.io/allen-discovery-center-mcovert/wcm-code:latest

RUN python runscripts/manual/runParca.py

CMD ["/bin/bash"]
