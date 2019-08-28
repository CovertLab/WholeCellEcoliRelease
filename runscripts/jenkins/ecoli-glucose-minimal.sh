HOST=$1
NAME=$2
PORT=$3
PASSWORD=$4

set -e

runscripts/jenkins/purge.sh daily_build 10

module load wcEcoli/sherlock2
pyenv local wcEcoli2

make clean
make compile

sh runscripts/jenkins/fireworks-config.sh $HOST $NAME $PORT $PASSWORD

echo y | lpad reset

PYTHONPATH=$PWD DESC="Daily build." SINGLE_DAUGHTERS=1 N_GENS=25 MASS_DISTRIBUTION=0 COMPRESS_OUTPUT=1 PLOTS=ACTIVE RAISE_ON_TIME_LIMIT=1 python runscripts/fireworks/fw_queue.py

PYTHONPATH=$PWD rlaunch rapidfire --nlaunches 0

N_FAILS=$(lpad get_fws -s FIZZLED -d count)

if [ $N_FAILS -gt 0 ]; then
  mv out/2* /scratch/PI/mcovert/wc_ecoli/failed/
fi

test $N_FAILS = 0

cp out/2*/kb/rawData.cPickle.bz2 /scratch/PI/mcovert/wc_ecoli/cached/
bunzip2 -f /scratch/PI/mcovert/wc_ecoli/cached/rawData.cPickle.bz2
chmod 444 /scratch/PI/mcovert/wc_ecoli/cached/rawData.cPickle
cp out/2*/kb/simData.cPickle.bz2 /scratch/PI/mcovert/wc_ecoli/cached/
bunzip2 -f /scratch/PI/mcovert/wc_ecoli/cached/simData.cPickle.bz2
chmod 444 /scratch/PI/mcovert/wc_ecoli/cached/simData.cPickle

mv out/2* /scratch/PI/mcovert/wc_ecoli/daily_build/
