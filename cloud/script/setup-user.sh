#!/usr/bin/env bash

set -eu

NAME=$1
EMAIL=$2

echo ">>>> establishing git configuration"
git config --global user.name "$NAME"
git config --global user.email "$EMAIL"

git config --global alias.co checkout
git config --global alias.br branch
git config --global alias.ci commit
git config --global alias.st status

echo ">>>> generating ssh key for $EMAIL"
ssh-keygen -t rsa -b 4096 -C "$EMAIL"
cat ~/.ssh/id_rsa.pub
