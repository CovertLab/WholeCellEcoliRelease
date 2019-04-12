# If you have docker installed, you can build this container with the following command:
#
#     > docker build -t wholecell .
#
# Then you can enter the container with the following command:
#
#     > docker run -it -d wholecell --name wholecelltest /bin/sh
#
# You will be dropped into an `sh` shell where you can execute commands:
#
#     # nosetests -a smalltest
#
# If this succeeds you should be good to go. Next, run the Parca:
#
#     # python runscripts/manual/runParca.py wholecellcontainer


FROM python:2.7.15-alpine

RUN apk add --no-cache build-base gcc make cmake g++ wget curl llvm pkgconfig freetype-dev libpng-dev swig openblas-dev gfortran ncurses-dev librdkafka-dev linux-headers

RUN wget ftp://ftp.gnu.org/gnu/glpk/glpk-4.65.tar.gz
RUN tar -xzvf glpk-4.65.tar.gz
RUN cd /glpk-4.65 && ./configure && make install && cd ..

COPY requirements.txt /
RUN pip install 'numpy==1.14.5'
RUN CVXOPT_BUILD_GLPK=1 pip install -r requirements.txt

COPY ./ /wcEcoli
WORKDIR /wcEcoli

RUN make clean compile
RUN export PYTHONPATH=`pwd`