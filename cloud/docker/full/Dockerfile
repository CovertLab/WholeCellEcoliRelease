# Container image #3: wcm-full.
# This Dockerfile builds a container image with the full Whole Cell Model, ready
# to run simulations, layered on the wcm-code image.
#
# To build locally from the wcEcoli/ project root directory:
#
#     > docker build -f cloud/docker/full/Dockerfile -t wcm-full .
#
# After building locally you can start up a new container from the image:
#
#     > docker run --name wcm -it --rm wcm-full
#
# It will start a shell where you can execute commands:
#
#     # nosetests
#
# If this succeeds you should be good to go, e.g.:
#
#     # python runscripts/manual/runSim.py


ARG from=wcm-code:latest
FROM ${from}

RUN python runscripts/manual/runFitter.py

CMD ["/bin/bash"]
