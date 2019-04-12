FROM python:2.7.15-alpine

RUN apk add --no-cache build-base gcc make cmake g++ wget curl llvm pkgconfig freetype-dev libpng-dev swig openblas-dev gfortran ncurses-dev librdkafka-dev linux-headers

RUN wget ftp://ftp.gnu.org/gnu/glpk/glpk-4.65.tar.gz
RUN tar -xzvf glpk-4.65.tar.gz
RUN cd /glpk-4.65 && ./configure && make install && cd ..

COPY requirements.txt /
RUN pip install 'numpy==1.14.5'
RUN pip install -r requirements.txt

COPY ./ /wcEcoli
WORKDIR /wcEcoli

RUN make clean compile
RUN export PYTHONPATH=`pwd`