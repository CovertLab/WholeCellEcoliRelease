# Container image #2: wcm-code.
# This Dockerfile builds a container image with the wcEcoli whole cell model
# code, layered on the wcm-runtime container image.
#
# To build this image locally from the wcEcoli/ project root directory:
#
#     > docker build -f cloud/docker/wholecell/Dockerfile -t wcm-code .
#
# (Add option `--build-arg from=ABC` to build from a different base image.)
#
# (If you want to build container image #3 on top of the locally-built base,
# you'll have to edit that Dockerfile or upload this image.)
#
# (To build using the Cloud Build service and store in the Container Registry,
# run `cloud/build.sh`.)
#
# After building locally you can start up a new container from the image:
#
#     > docker run --name wholecelltest -it wcm-code
#
# or if you used build.sh or build-wcm.sh to build using Cloud Build:
#
#     > docker run --name wholecelltest -it gcr.io/allen-discovery-center-mcovert/$USER-wcm-code
#
# It will start a shell where you can execute commands:
#
#     # nosetests
#
# If this succeeds you can start running WCM code in the container, e.g.:
#
#     # python runscripts/manual/runParca.py

ARG from=wcm-runtime:latest
FROM ${from}

LABEL application="Whole Cell Model of Escherichia coli" \
    email="allencentercovertlab@gmail.com" \
    license="https://github.com/CovertLab/WholeCellEcoliRelease/blob/master/LICENSE.md" \
    organization="Covert Lab at Stanford" \
    website="https://www.covert.stanford.edu/"

COPY . /wcEcoli
WORKDIR /wcEcoli

RUN make clean compile
ENV PYTHONPATH=/wcEcoli

RUN echo 'cat cloud/docker/welcome.txt' >> $HOME/.bashrc

CMD ["/bin/bash"]
