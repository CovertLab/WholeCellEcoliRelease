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
ARG from=python:3.8.3
FROM ${from}

# Add option `--build-arg NO_AVX2=1` to make OpenBLAS 0.3.6+ avoid AVX2
# instructions so it can pass its build self-test in Docker-for-Mac and
# produce correct results when run there.
#
# * NO_AVX2=1 is needed to build or run properly in Docker-for-Mac.
# * NO_AVX2=0 reportedly runs 20-30% faster BLAS if you're sure it'll only
#   build and run on Linux, but this only saves ~7% in a cell sim.
#
# Docker-for-Mac bug:
# https://github.com/xianyi/OpenBLAS/issues/2244
# https://github.com/docker/for-mac/issues/4576
# https://github.com/machyve/xhyve/issues/171
ARG NO_AVX2=1
ENV NO_AVX2="$NO_AVX2"

RUN apt-get update \
    && apt-get install -y swig gfortran llvm cmake nano

RUN echo "alias ls='ls --color=auto'" >> ~/.bashrc \
    && echo "alias ll='ls -l'" >> ~/.bashrc

# Install openblas (must be 0.3.5+, and 0.3.6+ need "NO_AVX2=1" for macOS Docker
# to produce correct results) in Theano and maybe numpy/scipy.
# Someday `apt-get install libopenblas-dev=...` will work.
ENV OPENBLAS_LABEL=v0.3.10
RUN (mkdir -p openblas && cd openblas \
    && curl -SL https://github.com/xianyi/OpenBLAS/archive/${OPENBLAS_LABEL}.tar.gz | tar -xz \
    && cd OpenBLAS* \
    && echo "Compiling OpenBLAS ${OPENBLAS_LABEL} with NO_AVX2=${NO_AVX2}" \
    && make "NO_AVX2=${NO_AVX2}" FC=gfortran \
    && make "NO_AVX2=${NO_AVX2}" PREFIX=/usr install) \
    && rm -r openblas

# Conditionally build numpy and scipy from source (`--no-binary`, below) to
# link them to this OpenBLAS instead of their embedded copies of OpenBLAS.
# That can add 50 minutes to the build time.
RUN if [ "$NO_AVX2" != 0 ] ; then (echo "[openblas]" \
    && echo "libraries = openblas" \
    && echo "library_dirs = /usr/lib" \
    && echo "include_dirs = /usr/include") > ~/.numpy-site.cfg ; \
    fi

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
    && pip install --no-cache-dir numpy==1.19.0 $b1 \
    && pip install --no-cache-dir -r requirements.txt $b1 \
    && umask 000 && mkdir -p /.theano)

CMD ["/bin/bash"]
