# This is the Dockerfile for building the environment for the whole cell model. It stops just
# short of copying in the `wcEcoli` code.
#
# To build locally, execute the following command:
#
#     > docker build -t wcm-runtime .
#
# To build using gcloud container service:
#
#     > gcloud builds submit --timeout=2h --tag gcr.io/allen-discovery-center-mcovert/wcm-runtime .

FROM python:2.7.15-alpine

RUN apk add --no-cache build-base linux-headers gcc g++ gfortran llvm make cmake wget curl pkgconfig freetype-dev libpng-dev swig openblas-dev ncurses-dev librdkafka-dev

RUN wget ftp://ftp.gnu.org/gnu/glpk/glpk-4.65.tar.gz
RUN tar -xzvf glpk-4.65.tar.gz
RUN cd /glpk-4.65 && ./configure && make install && cd ..

ENV CVXOPT_BUILD_GLPK=1
ENV CVXOPT_GLPK_LIB_DIR=/usr/local/lib
ENV CVXOPT_GLPK_INC_DIR=/usr/local/include
ENV CVXOPT_BLAS_LIB=openblas
ENV CVXOPT_LAPACK_LIB=openblas

# Install all the pips within one Docker layer and don't cache the downloads.
COPY requirements.txt /
RUN pip install --no-cache-dir 'numpy==1.14.5' \
    && pip install --no-cache-dir -r requirements.txt
