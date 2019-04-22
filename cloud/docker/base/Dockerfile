# This is the Dockerfile for building the environment for the whole cell model. It stops just
# short of copying in the `wcEcoli` code.
#
# To build locally, execute the following command in the root level wcEcoli directory:
#
#     > docker build -f cloud/docker/base/Dockerfile -t wcm-runtime .
#
# To build using gcloud container service in this directory:
#
#     > gcloud builds submit --timeout=2h --tag gcr.io/allen-discovery-center-mcovert/wcm-runtime .

FROM python:2.7.15-alpine

RUN apk add --no-cache build-base linux-headers gcc g++ gfortran llvm make cmake \
    wget curl pkgconfig freetype-dev libpng-dev swig ncurses-dev librdkafka-dev

# Install openblas v0.3.5 (not buggy v0.3.3).
# Future:  RUN apk add --no-cache openblas-dev
RUN apk add --no-cache perl
RUN (mkdir -p openblas && cd openblas \
    && curl -SL https://github.com/xianyi/OpenBLAS/archive/v0.3.5.tar.gz | tar -xz \
    && cd OpenBLAS* \
    && make FC=gfortran \
    && make PREFIX=/usr install) \
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
    CVXOPT_LAPACK_LIB=openblas

# Install all the pips within one Docker layer and don't cache the downloads.
COPY requirements.txt /
RUN pip install --no-cache-dir 'numpy==1.14.5' \
    && pip install --no-cache-dir -r requirements.txt
