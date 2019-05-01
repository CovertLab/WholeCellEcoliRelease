#!/usr/bin/env bash

set -eu

DIR=$1
DEVICE=$2

echo ">>>> formatting the disk /dev/$DEVICE"
mkfs.ext4 -m 0 -F -E lazy_itable_init=0,lazy_journal_init=0,discard /dev/$DEVICE

echo ">>>> making mount directory at /mnt/disks/$DIR"
mkdir -p /mnt/disks/$DIR

echo ">>>> mounting /dev/$DEVICE at /mnt/disks/$DIR"
mount -o discard,defaults /dev/$DEVICE /mnt/disks/$DIR

echo ">>>> making /mnt/disks/$DIR writable"
chmod a+w /mnt/disks/$DIR
