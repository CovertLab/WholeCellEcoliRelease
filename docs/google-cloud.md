# How to run the Whole Cell Model on the Google Cloud Platform

## Setup

* One time only: [Install and log in to the Google Cloud SDK](https://github.com/CovertLab/sisyphus/blob/master/GCLOUD_SETUP.md).
* The steps below assume your current working directory is your wcEcoli git clone:  
  `cd ~/dev/wcEcoli`
* The `python` steps assume the wcEcoli directory is on the Python path:  
  `export PYTHONPATH=$PWD`
* In your web browser, bookmark the [Google Cloud Platform
console](https://console.cloud.google.com/home/dashboard?project=allen-discovery-center-mcovert)


## Run it

1. Upload and build a Docker image containing the Whole Cell Model code:

   ```sh
   cloud/build-wcm.sh
   ```

   Notes:

   * This creates a Docker image named `gcr.io/allen-discovery-center-mcovert/$USER-wcm-code`,
   where `$USER` is your local username.
   * It takes 1/2 minute to upload the source code and 2 minutes for the Google Cloud Build
   server to build the Docker container image including compiling the project's Cython code.
   * **Time saver:** _You can skip this step_ if the Whole Cell Model code (Parca, simulation, and
   analysis) didn't change since you last did it.
   * See "Variations" below about multiple `wcm-code` Docker images and new `wcm-runtime` base images.

2. Open an ssh connection and TCP tunnel to the Gaia workflow server:

   ```sh
   runscripts/cloud/ssh-tunnel.sh
   ```

3. Use another terminal tab to build, upload, and start the Whole Cell Model workflow:

   ```sh
   python runscripts/cloud/wcm.py
   ```

   Notes:

   * Use the `-h` option to get help on all the command line options, e.g. multiple
   simulation generations.
   * Use the `--dump` option to write the workflow definition as JSON files for
   review instead of uploading them to the server. You can edit these files and
   upload them manually or re-run `wcm.py` without `--dump`.
   * After `wcm.py` runs, you can quit the `ssh-tunnel.sh` unless you want to
   use it for interactive Gaia operations in a Python terminal.


## Watch it run

[TODO]


## Download the outputs

[TODO]


## Variations

### A new `wcm-runtime` base Docker image

[TODO]

Times for building wcm-runtime: It takes 35:40 minutes:seconds of installing stuff plus docker I/O work. Out of that, 0:43 for apt-get, 13:56 to install openblas, 3:09 to install numpy, and 14:51 to install scipy (after downloading it along with the other pips)


### Multiple `wcm-code` Docker images

If you want to build an independent `wcm-code` Docker image, pass an optional "workflow ID" to
`cloud/build-wcm.sh`, e.g.:

   ```sh
   cloud/build-wcm.sh nightly
   ```

[TODO] document how to pass this workflow ID to `wcm.py`.


### Uploading workflow definition files manually

[TODO]


### Future

* Test other ways to open a secure connection to the Gaia workflow server.
