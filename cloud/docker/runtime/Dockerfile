# Container image #1: wcm-runtime.
# This Dockerfile builds the runtime environment for the whole cell model.
#
# To build this image locally from the wcEcoli/ project root directory:
#
#     > docker build -f cloud/docker/runtime/Dockerfile -t ${USER}-wcm-runtime .
#
# (To build using the Cloud Build service and store in the Container Registry,
# run `cloud/build.sh`.)
#
# Add option `--build-arg from=ABC` to build from a different base image "ABC"
# but DO NOT USE an alpine base since the simulation math comes out different!
# See https://pythonspeed.com/articles/alpine-docker-python/ for more reasons.
ARG from=python:3.8.5
FROM ${from}

RUN echo "alias ls='ls --color=auto'" >> ~/.bashrc \
    && echo "alias ll='ls -l'" >> ~/.bashrc \
    && cp ~/.bashrc /

# Update and install in the same layer so it won't install from old updates.
RUN apt-get update \
    && apt-get install -y swig gfortran llvm cmake nano

# Option `--build-arg NO_AVX2=1` controls OpenBLAS' use of AVX2 vector instructions:
#
# * NO_AVX2=1 is needed to build properly in Docker-for-Mac due to a Docker bug,
#   but it computes slightly different results.
# * NO_AVX2=0 gets more consistent results and reportedly runs 20-30% faster
#   (although it only saves ~7% in a cell sim). Use this when building in
#   Docker-for-Linux. The built Images seem to run fine on Mac.
#
# Docker-for-Mac bug:
# https://github.com/xianyi/OpenBLAS/issues/2244
# https://github.com/docker/for-mac/issues/4576
# https://github.com/machyve/xhyve/issues/171
ARG NO_AVX2=0
ENV NO_AVX2="$NO_AVX2"

# Install openblas (must be 0.3.5+) for numpy, scipy, and Theano.
ENV OPENBLAS_LABEL=v0.3.10
RUN (mkdir -p openblas && cd openblas \
    && curl -SL https://github.com/xianyi/OpenBLAS/archive/${OPENBLAS_LABEL}.tar.gz | tar -xz \
    && cd OpenBLAS* \
    && echo "Compiling OpenBLAS ${OPENBLAS_LABEL} with NO_AVX2=${NO_AVX2}" \
    && make "NO_AVX2=${NO_AVX2}" FC=gfortran \
    && make "NO_AVX2=${NO_AVX2}" PREFIX=/usr install) \
    && rm -r openblas

# Build numpy and scipy from source (`--no-binary`) to link them to this
# OpenBLAS. This can add 50 minutes to the build time to get more consistent
# results and go faster at run time.
RUN (echo "[openblas]" \
    && echo "libraries = openblas" \
    && echo "library_dirs = /usr/lib" \
    && echo "include_dirs = /usr/include") > ~/.numpy-site.cfg

# This gets more consistent results from openblas.
ENV OPENBLAS_NUM_THREADS=1

# Install all the pips within one Docker layer, don't cache the downloads, and
# enable Theano to write into the data dir it gets when running as a user with
# no home dir. (That happens because DockerTask runs tasks inside the Docker
# container with the uid:gid of the user outside the container to get
# read/write/delete access to the tasks's output files.)
COPY requirements.txt /
RUN (b1="" \
    && if [ -f ~/.numpy-site.cfg ] ; then b1="--no-binary=numpy,scipy"; fi \
    && echo "Installing pips with $b1" \
    && pip install --no-cache-dir --upgrade pip setuptools wheel \
    && pip install --no-cache-dir numpy==1.19.2 $b1 \
    && pip install --no-cache-dir -r requirements.txt $b1 \
    && umask 000 && mkdir -p /.theano)

CMD ["/bin/bash"]
