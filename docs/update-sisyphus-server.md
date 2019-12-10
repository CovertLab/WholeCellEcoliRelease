# How to update the Sisyphus server's disk Image

Sisyphus servers on Google Compute Engine (GCE) are the "worker bees" that run
the workflow "payload" steps.

It's good to occasionally update the disk image that they start from to get
security updates and to preload the current wcm-runtime base Docker image,
which saves a bunch of seconds every time a worker instance starts up.


## Setup

This assumes you installed the `gcloud` command line tool then logged in and
configured it via `gcloud init`. Check the settings:

    gcloud config list

the printout should include:

        [compute]
        region = us-west1
        zone = us-west1-b
        [core]
        account = «YOUR LOGIN ACCOUNT»
        project = allen-discovery-center-mcovert

If those settings aren't right, you can go through `gcloud init` again or quickly set

    gcloud config set compute/zone us-west1-b
    gcloud config set core/project allen-discovery-center-mcovert


## Steps

* If you need to increase the size of the Sisyphus boot disk, edit the
`--boot-disk-size=200GB` line in `runscripts/cloud/launch-workers.sh`.

* Launch an instance of the server

      WORKFLOW=none runscripts/cloud/launch-workers.sh resizing

  (We normally launch these servers with a name that starts with "sisyphus-".
  Here we use an unusual name to find it easily in the Google Cloud Console
  and let people filter it out in the Logs Viewer.)

* ssh to it and stop the sisyphus worker service so it won't change things, or
start running workflow steps, or delete itself

      gcloud compute ssh resizing
      sudo systemctl stop sisyphus.service

* If the OS printed `*** System restart required ***` on login, do

      sudo reboot

  then `ssh` to it again and `sudo systemctl stop` the service again.

* If you're resizing the boot disk, verify that it has the new disk size via
`sudo lsblk` and `df -h .` (This OS should auto-resize the partition table, so
you don't have to mess with that.)

* Update the apt packages

      sudo apt-get update
      sudo apt-get upgrade  # if "Resource temporarily unavailable", wait and try again
      sudo apt autoremove
      sudo reboot

  then `ssh` to it again and `sudo systemctl stop` the service again.

* Trim the service journal so we don't keep seeing old failures and junk when
using `sudo journalctl -u sisyphus` or `sudo journalctl -u sisyphus -f` on the
server to view the service logs

      sudo journalctl --vacuum-time=2d


## Do the following steps as the sisyphus user

      sudo su -l sisyphus

* Clean out all Docker containers and old (or all) docker images, if any

      docker ps -aq
      docker stop $(docker ps -aq)
      docker rm $(docker ps -aq)

* Delete old Docker images [maybe keep the immediately previous `wcm-runtime` image]

      docker image ls --all
      docker image rm gcr.io/allen-discovery-center-mcovert/wcm-runtime

* Pull the latest Docker images `wcm-runtime` and `python:2.7.16`

      docker pull python:2.7.16
      docker pull gcr.io/allen-discovery-center-mcovert/wcm-runtime

* Pull the latest Sisyphus git source code from master

      cd
      cd sisyphus
      git pull origin master

* Exit from the sisyphus user account

      ^D


## Now make the disk image

* Stop the server so its Persistent Disk is in a clean state to make an Image

      sudo shutdown -h now

* Wait for it to be stopped, checking the Compute Engine > VM instances page

* Look in the Compute Engine > Images page for the last disk image name in the `sisyphus-v*` series

* In the Compute Engine > Disks page
  1. click on the Disk named `resizing`
  2. click `CREATE IMAGE`
  3. fill these properties into the form:

         Name: sisyphus-v6             # whatever's the next version number in the series
         Source: Disk
         Source disk: resizing
         Location: Multi-regional      # [change it to Regional???]
         Family: sisyphus-worker       # <====== MUST DO ======
         Description: Sisyphus Worker  # also list the Docker images and sisyphus git hash
         Labels: role=worker

  4. click the `Create` button to make a new disk Image
  5. watch the Console to see when it's done

* Test it by running the "demo" workflow in the usual way

      python runscripts/cloud/util/demo_workflow.py

* In the Console, delete old Images in this Family, keeping the one or two previous Images

* Delete the `resizing` GCE Instance (and its boot disk, too)


## Notes

* If the new Disk Image doesn't work, you can mark it "deprecated" in the Image Family to
revert to the previous Image.
* Rebooting a GCE Instance will fix the `*** System restart required ***` login message but
it won't fix the message about _n_ packages that need updating.
* Restarting a server seems to preserve its internal IP address and change its external IP
address. For local connections between servers in the data center, use the instance names
which are the same as their boot Disk names: `gaia-base`, `rabbit-prime`, and `zookeeper-prime`.
