HOST=$1
NAME=$2
PORT=$3
PASSWORD=$4

set -e

mkdir -p /scratch/PI/mcovert/jenkins/fireworks/logs/launchpad
echo "logdir: /scratch/PI/mcovert/jenkins/fireworks/logs/launchpad" >> my_launchpad.yaml
echo "host: $HOST" >> my_launchpad.yaml
echo "name: $NAME" >> my_launchpad.yaml
echo "username: fireworks" >> my_launchpad.yaml
echo "password: $PASSWORD" >> my_launchpad.yaml
echo "port: $PORT" >> my_launchpad.yaml
echo "strm_lvl: INFO" >> my_launchpad.yaml
echo "user_indices: []" >> my_launchpad.yaml
echo "wf_user_indices: []" >> my_launchpad.yaml

# So pymongo 3.9.0+ will work with older MongoDB servers:
echo "mongoclient_kwargs:" >> my_launchpad.yaml
echo "  retryWrites: false" >> my_launchpad.yaml
