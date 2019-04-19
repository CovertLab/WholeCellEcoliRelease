# If you have docker installed, you can build this image with the following command:
#
#     > docker build -t wholecell -f cloud/docker/wholecell/Dockerfile .
#
# Then you can create a new container from this image with the following command:
#
#     > docker run --name wholecelltest -it -d wholecell
#
# Once the container exists you can `exec` it to get a shell inside:
#
#     > docker exec -it wholecelltest /bin/sh
#
# You will be dropped into an `sh` shell where you can execute commands:
#
#     # nosetests
#
# If this succeeds you should be good to go. Next, run the Parca:
#
#     # python runscripts/manual/runParca.py


FROM gcr.io/allen-discovery-center-mcovert/wcm-runtime:latest

COPY . /wcEcoli
WORKDIR /wcEcoli

RUN make clean compile
ENV PYTHONPATH=/wcEcoli
