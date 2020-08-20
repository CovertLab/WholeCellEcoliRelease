# Container image #1: wcm-runtime.
# This Dockerfile builds the runtime environment for the whole cell model.
#
# To build this image locally from the wcEcoli/ project root directory:
#
#     > docker build -f cloud/docker/runtime/Dockerfile -t wcm-runtime .
#
# After building locally you can start up a new container from the image:
#
#     > docker run --name runtime-test -it --rm wcm-runtime
#
# It will start a shell where you can execute commands:
#
#     # python -c 'import numpy; numpy.show_config()'
#     # python -c 'import scipy; scipy.__config__.show()'
#
# Add option `--build-arg from=ABC` to build from a different base image but
# DO NOT USE an alpine base since the simulation math comes out different!
# See https://pythonspeed.com/articles/alpine-docker-python/ for more reasons.
ARG from=python:2.7.16-stretch
FROM ${from}

# Option `--build-arg NO_AVX2=1` controls OpenBLAS' use of AVX2 instructions:
#
# * NO_AVX2=1 is needed to build properly in Docker-for-Mac due to a
#   Docker bug.
# * NO_AVX2=0 reportedly runs 20-30% faster BLAS if you're sure it'll only
#   build in Docker-for-Linux, but this only saves ~7% in a cell sim.
#
# Docker-for-Mac bug:
# https://github.com/xianyi/OpenBLAS/issues/2244
# https://github.com/docker/for-mac/issues/4576
# https://github.com/machyve/xhyve/issues/171
ARG NO_AVX2=0
ENV NO_AVX2="$NO_AVX2"

RUN apt-get update \
    && apt-get install -y swig gfortran llvm cmake nano python-cvxopt

RUN echo "alias ls='ls --color=auto'" >> ~/.bashrc \
    && echo "alias ll='ls -l'" >> ~/.bashrc

# Install openblas.
ENV OPENBLAS_LABEL=v0.3.9
RUN (mkdir -p openblas && cd openblas \
    && curl -SL https://github.com/xianyi/OpenBLAS/archive/${OPENBLAS_LABEL}.tar.gz | tar -xz \
    && cd OpenBLAS* \
    && echo "Compiling OpenBLAS ${OPENBLAS_LABEL} with NO_AVX2=${NO_AVX2}" \
    && make "NO_AVX2=${NO_AVX2}" FC=gfortran \
    && make "NO_AVX2=${NO_AVX2}" PREFIX=/usr install) \
    && rm -r openblas

# Install glpk into /usr/local/lib/.
RUN (mkdir -p glpk && cd glpk \
    && curl -SL ftp://ftp.gnu.org/gnu/glpk/glpk-4.65.tar.gz | tar -xz \
    && cd glpk-4.65 \
    && ./configure && make install) \
    && rm -r glpk

ENV CVXOPT_BUILD_GLPK=1 \
    CVXOPT_GLPK_LIB_DIR=/usr/local/lib \
    CVXOPT_GLPK_INC_DIR=/usr/local/include \
    CVXOPT_BLAS_LIB=openblas \
    CVXOPT_LAPACK_LIB=openblas \
    OPENBLAS_NUM_THREADS=1

# Install all the pips within one Docker layer and don't cache the downloads, and
# enable Theano to write into the data dir it gets when running as a user with
# no home dir. (That happens when running the process inside the Docker container
# with the uid:gid of the user outside the container in order to provide
# read/write/delete access to the process's output files.)
COPY requirements.txt /
RUN pip install --no-cache-dir 'numpy==1.14.6' --no-binary numpy \
    && pip install --no-cache-dir -r requirements.txt --no-binary numpy,scipy,cvxopt \
    && umask 000 && mkdir -p /.theano

CMD ["/bin/bash"]
