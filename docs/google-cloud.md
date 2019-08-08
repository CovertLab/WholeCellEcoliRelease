# How to run the Whole Cell Model on the Google Cloud Platform

## Setup

* One time setup: [Install and log in to the Google Cloud SDK](https://github.com/CovertLab/sisyphus/blob/master/GCLOUD_SETUP.md).
* The steps below assume your current working directory is your wcEcoli git clone, e.g.:  
  `cd ~/dev/wcEcoli`
* The Python steps assume the wcEcoli directory is on the Python path:  
  `export PYTHONPATH=$PWD`


## Run it

1. Upload your Whole Cell Model code:

   ```sh
   cloud/build-wcm.sh
   ```

   Notes:

   * This uploads the source code from your working directory, whether or not it's checked-in or
   compiled, so consider running small-scale local tests beforehand using
   [the manual runscripts](https://github.com/CovertLab/wcEcoli#using-the-manual-runscripts).
   You can shorten the run time via command line options like runSim's `--length-sec`.
   * This creates a Docker image named `gcr.io/allen-discovery-center-mcovert/$USER-wcm-code` in
   the project's section of Google Cloud Registry, where `$USER` is your local username.
   * It takes 1/2 minute to upload the source code and 2 minutes for the Google Cloud Build
   server to build the Docker container image, including compiling the project's Cython code.
   * **Time saver:** _You can skip this step_ if the Whole Cell Model code (Parca, simulation, and
   analysis) didn't change since you last ran this step.
   * See "Variations" below about multiple `wcm-code` Docker images and new `wcm-runtime` base images.

2. Open an ssh connection and TCP tunnel to the Gaia workflow server:

   ```sh
   runscripts/cloud/ssh-tunnel.sh
   ```

3. Use another terminal tab to build and start the Whole Cell Model workflow:

   ```sh
   python runscripts/cloud/wcm.py
   ```

   Notes:

   * Use the `-h` option to get help on all the command line options, e.g. how many
   simulation generations to run, how many initial cell seeds, and a workflow
   description.
   * Use the `--dump` option to write the workflow definition as JSON files for
   review instead of running them in Google Compute Engine. You can edit these
   files and manually upload them to the Gaia workflow server or just re-run
   `wcm.py` without `--dump`.
   * After `wcm.py` runs, you can quit the `ssh-tunnel.sh` unless you want to
   use it for interactive Gaia operations in a Python terminal.
   * `wcm.py` will print the workflow's timestamped name, e.g.
   `WCM_$USER_20190725.160911`, with your local user name.
   Watch for this workflow name in the logs. 


## Watch it run

### Setup

* In your web browser, bookmark the [Google Cloud Platform
console](https://console.cloud.google.com/home/dashboard?project=allen-discovery-center-mcovert).
* In the `☰` menu (pronounced "hamburger") in the top left corner of the Google Cloud Platform console webpage,
pin "Compute Engine" (GCE), "Logging" (Stackdriver), and "Storage" (GCS) to the top section of this menu
for quick access.
* You can rearrange and hide the "cards" on this page.

**Optional:** Set up a chart to display worker node CPU utilization.
(_"Worker nodes"_ are the Compute Engine virtual machines that run the "real
work" steps of the simulation and analysis.)

  1. Move the "Compute Engine" card to top dead center position.
  2. In this card's `⋮` menu, click `Add Chart`.
  3. Title it `Worker load`.
  4. Add Metric type `instance/cpu/utilization` with
  filter `metric.labels.instance_name = starts_with("sisyphus")`.
  5. Add Metric type `instance/disk/read_bytes_count` with the same filter,
  group by `resource.label.project_id`, aggregate `Sum`.
  5. Add Metric type `instance/disk/write_bytes_count` with the same filter,
  group by `resource.label.project_id`, aggregate `Sum`.
  6. Click `Save`.



### Monitor

* Click your web browser bookmark to open the [Google Cloud Platform
console](https://console.cloud.google.com/home/dashboard?project=allen-discovery-center-mcovert)
home page. The `☰` menu is a way to navigate to all the other pages.

* Open the [Compute Engine — VM instances](https://console.cloud.google.com/compute/instances?project=allen-discovery-center-mcovert&instancessize=50)
page to see the list of running Compute Engine VM instances.

   * "sisyphus-…" VM instances are the "worker bees" that run the workflow steps; `wcm.py`
   launches workers named `sisyphus-$USER-0`, `sisyphus-$USER-1`, … but (for now) they are
   not dedicated to your workflow
   * "gaia-base" is always available to manage workflows
   * "rabbit-prime" is running a queue of workflow tasks
   * "zookeeper-prime" communicates control information between Gaia and worker nodes

* Open the [Logging — Logs Viewer](https://console.cloud.google.com/logs/viewer?project=allen-discovery-center-mcovert)
page to view the logs from the project's GCE VM instances.

   * The triangle "Run" button at the top starts streaming the logs (akin to `tail -f` following).
   * The page has several filtering tools. A good place to start is to set the resource menu to
   `GCE VM Instance` and the log level to `Info`.
   * Log level `Debug` will show internal workings of the Gaia and Sisyphus servers.
   * Log severity `NOTICE` is used to make startup and shutdown events more prominent
   than `INFO`, but Log Viewer shows them with the same icon as `INFO` entries.
   * You can also filter on your workflow name or more simply on your user name.
   * Log entries from the Gaia server have
   logName = `projects/allen-discovery-center-mcovert/logs/gaia-base` and
   resource type `gce_instance` and resource label instance_id = `gaia-base`.
   * Log entries from Sisyphus servers have logNames like
   `projects/allen-discovery-center-mcovert/logs/sisyphus-$USER-0.Demo_$USER_20190729.193458.count`
   when running a step (named `count` in this example) of a workflow (named `Demo_$USER_20190729.193458`)
   in this example, or just `projects/allen-discovery-center-mcovert/logs/sisyphus-$USER-0`
   when not running a workflow step. They have
   resource type `gce_instance` and resource label instance_id = `sisyphus-$USER-$NUMBER`.
   * The logs should show all the output lines from your workflow steps and some log entries
   from Gaia and Sisyphus coordination.
   * When it logs `WORKFLOW COMPLETE` naming your workflow, your workflow is done!

* Open the [Storage — Browser](https://console.cloud.google.com/storage/browser?project=allen-discovery-center-mcovert)
page to browse the files created by our workflow runs.

   * Currently the files are all within Bucket `sisyphus`, directory `data`.
   * `wcm.py` directs its output files ito a subdirectory named for the user, the
   date-time timestamp, and the optional workflow description.


## Download the outputs

There are many ways to download the outputs from your workflow:

* **Simplest:** Open the
[Storage — Browser](https://console.cloud.google.com/storage/browser?project=allen-discovery-center-mcovert),
browse into the [sisyphus storage bucket](https://console.cloud.google.com/storage/browser/sisyphus?project=allen-discovery-center-mcovert)
bucket, find your workflow files, and click on individual files to download them.
* **Most convenient:** Use [gcsfuse](https://github.com/GoogleCloudPlatform/gcsfuse) to mount the storage
bucket `sisyphus` to your local file system (or mount just its `data/$USER/` subdirectory) and find
the files you want.
  * gcsfuse reads and writes whole files to Cloud Storage on demand. It's convenient but
  it has high latency.
  * Cloud Storage does not act like a regular file system and gcsfuse cannot hide that.
  E.g. it doesn't have directories, just file paths that may contain slashes. See the notes in
  [semantics.md](https://github.com/GoogleCloudPlatform/gcsfuse/blob/master/docs/semantics.md)
  for details.
  * For now, you'll need to use gcsfuse's `--implicit-dirs` option until Sisyphus creates stand-ins
  for directory entries. This makes gcsfuse even slower and has
  [additional side-effects](https://github.com/GoogleCloudPlatform/gcsfuse/blob/master/docs/semantics.md#implicit-directories).
* **Fastest:** Use the [gsutil tool](https://cloud.google.com/storage/docs/gsutil) (it's part of the
Cloud SDK tools you installed along with `gcloud`) to do anything you need with Cloud Storage
buckets and files. With the `-m` option, it'll transfer lots of files in parallel, so this is the
fastest approach.
  * It can also do [gsutil -m rsync](https://cloud.google.com/storage/docs/gsutil/commands/rsync).
* We might write a visual browser for pictures in our Cloud Storage bucket.
* If you have lots of microscopy images, try running [ffmpeg](https://ffmpeg.org/) to
[convert an image sequence to](https://trac.ffmpeg.org/wiki/Encode/AV1) an
[AV1](https://en.wikipedia.org/wiki/AV1) video file for much faster download and viewing. AV1
is a state of the art (released March 2018) video format that competes with H.265. It aims to
produce substantially better quality/size than H.264 and VP9 and also be royalty-free.
It does support lossless compression. It's supported by Chrome, Firefox, and Opera.


## Variations

### A new `wcm-runtime` base Docker image

If the WCM pips or apt-get packages change, use the `cloud/build-runtime.sh` script
to rebuild the `wcm-runtime` base Docker image, then build new `wcm-code` images
from that.

By default, we're all sharing a `wcm-runtime` base Docker image, so either ensure updates are
backwards compatible or pass a different WCM_RUNTIME name into the `cloud/build-wcm.sh` script.

Each Docker image layer has long hash references to its lower layer images (like a git commit
referring to the previous commit hash) so building a new `wcm-runtime` image won't alter existing
`wcm-code` images (unless someone gets overzealous about cleaning out images from the container registry).

When a new branch needs a new `wcm-runtime` image, pick a name for the new runtime like
`wcm-runtime:pr1234` (or maybe `wcm-runtime-pr1234`) and pass the name to the scripts:


   ```sh
   cloud/build-runtime.sh wcm-runtime:pr1234
   cloud/build-wcm.sh $USER wcm-runtime:pr1234
   python runscripts/cloud/wcm.py
   ```

or setting a workflow user ID to do a CI PR build:

   ```sh
   cloud/build-runtime.sh wcm-runtime:pr1234
   cloud/build-wcm.sh PR1234 wcm-runtime:pr1234
   WF_ID=PR1234 python runscripts/cloud/wcm.py
   ```

Q. After merging that PR into master, should we then build a new `wcm-runtime:latest` image,
changing the meaning of that default name for future `wcm-code` image builds? Or enshrine a
`:version` ID into the scripts?

Optimization: The Sisyphus Compute Engine disk image has a copy of the `wcm-runtime` Docker
image to save tens of seconds in startup time. After creating a new `wcm-runtime` Docker image,
it's a good idea to add that to the Sisyphus disk image.

Potential optimization: `cloud/build-runtime.sh` takes 36m installing stuff plus Docker I/O work.
Of that, 14m goes into installing openblas, 3m to installing numpy, and 15m to installing scipy
(after downloading it along with the other pips). So if we run this frequently enough, it would
be worth splitting out a first base image containing the apt-get steps, openblas, numpy, and
scipy from the `wcm-runtime` base image that adds all the other pips.

Beware that OpenBLAS changes seem to mostly be about compatibility with various
OS and other platform changes, and they're struggling at that. So OpenBLAS 0.3.5 works well
(although it's not available via apt-get) whereas 0.3.4 and 0.3.6 are broken.


### Multiple `wcm-code` Docker images

If you want to build more than one `wcm-code` Docker image so you can run different code in
different workflows, pass an optional "workflow user ID" to `cloud/build-wcm.sh`, e.g.:

   ```sh
   cloud/build-wcm.sh nightly
   ```

then pass the same workflow user ID to `wcm.py` via the environment variable `$WF_ID`. These
default to `$USER`.

   ```sh
   WF_ID=nightly python runscripts/cloud/wcm.py
   ```


### Future

* Document how to set up and update the Gaia, Sisyphus, RabbitMQ, and Kafka
servers on GCE.
  * Periodically security updates and apt-get upgrades.
  * Occasionally pull new Docker images to the Sisyphus server disk image and
  delete old ones.
  * Occasionally trim journalctl.
  * Sisyphus disk images are in an "image family" so new ones supersede older
  ones and we can revert back if a new one doesn't work. Do the same with the
  other server types. 
* Test other ways to open a secure connection to the Gaia workflow server.
* Make Sisyphus create the psuedo "directory" entries in Cloud Storage so we
can use `gcsfuse` without `--implicit-dirs`.
* A web UI for the Gaia workflow server to monitor and modify workflow runs.
