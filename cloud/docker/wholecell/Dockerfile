# Container image #2: wcm-code.
# This Dockerfile builds a container image with the wcEcoli whole cell model
# code, layered on the wcm-runtime container image.
#
# To build locally from the wcEcoli/ project root directory:
#
#     > docker build -f cloud/docker/wholecell/Dockerfile -t wcm-code .
#
# (If you want to build container image #3 on top of the locally-built base,
# you'll have to edit that Dockerfile or upload this image.)
#
# To build using the Cloud Build service and store in the Container Registry,
# from the wcEcoli/ project root directory:
#
#     > gcloud builds submit --timeout=15m --config config-build-2-wcm-code.json
#
# After building locally you can start up a new container from the image:
#
#     > docker run --name wholecelltest -it wcm-code
#
# It will start a shell where you can execute commands:
#
#     # nosetests
#
# If this succeeds you should be good to go, e.g.:
#
#     # python runscripts/manual/runParca.py


FROM gcr.io/allen-discovery-center-mcovert/wcm-runtime:latest

COPY . /wcEcoli
WORKDIR /wcEcoli

RUN make clean compile
ENV PYTHONPATH=/wcEcoli

CMD ["/bin/bash"]
