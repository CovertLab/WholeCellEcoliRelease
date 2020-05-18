# Container image #1: wcm-runtime.
# This Dockerfile builds the runtime environment for the whole cell model.
#
# To build this image locally from the wcEcoli/ project root directory:
#
#     > docker build -f cloud/docker/runtime/Dockerfile -t wcm-runtime .
#
# (To build using the Cloud Build service and store in the Container Registry,
# run `cloud/build.sh`.)
#
# Add option `--build-arg from=ABC` to build from a different base image but
# DO NOT USE an alpine base since the simulation math comes out different!
# See https://pythonspeed.com/articles/alpine-docker-python/ for more reasons.
ARG from=python:2.7.16
FROM ${from}

# Add option `--build-arg NO_AVX2=1` to make OpenBLAS 0.3.6+ not use the AVX2
# hardware so the image can build and run in Docker Desktop on macOS without bad
# results and fatal self-test errors.
# Turn this off once a Docker release has an xhyve release that fixes the bug.
# https://github.com/xianyi/OpenBLAS/issues/2244
# https://github.com/machyve/xhyve/issues/171
ARG NO_AVX2=0

RUN apt-get update \
    && apt-get install -y swig gfortran llvm cmake ncurses-dev \
        libreadline7 libreadline-dev nano

RUN echo "alias ls='ls --color=auto'" >> ~/.bashrc \
    && echo "alias ll='ls -l'" >> ~/.bashrc

# Install openblas (0.3.5 is OK; 0.3.3 is not; 0.3.6+ need "NO_AVX2=1"
# in macOS Docker to produce correct results).
# Someday `apt-get install libopenblas-dev=...` will work.
RUN (mkdir -p openblas && cd openblas \
    && curl -SL https://github.com/xianyi/OpenBLAS/archive/v0.3.9.tar.gz | tar -xz \
    && cd OpenBLAS* \
    && make "NO_AVX2=${NO_AVX2}" FC=gfortran \
    && make "NO_AVX2=${NO_AVX2}" PREFIX=/usr install) \
    && rm -r openblas

ENV OPENBLAS_NUM_THREADS=1

# Install all the pips within one Docker layer, don't cache the downloads, link
# numpy and scipy to the installed OpenBLAS, and enable Theano to write into
# the data dir it uses when running as a user with no home dir.
COPY requirements.txt /
RUN (echo "[openblas]" \
    && echo "libraries = openblas" \
    && echo "library_dirs = /usr/lib" \
    && echo "include_dirs = /usr/include") > ~/.numpy-site.cfg
RUN (pip install --no-cache-dir --upgrade pip setuptools wheel \
    && pip install --no-cache-dir numpy==1.14.6 --no-binary numpy \
    && pip install --no-cache-dir -r requirements.txt --no-binary numpy,scipy \
    && umask 000 && mkdir -p /.theano)

CMD ["/bin/bash"]
