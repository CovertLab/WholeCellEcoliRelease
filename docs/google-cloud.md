# How to run the Whole Cell Model on the Google Cloud Platform

➽ See [Borealis](https://github.com/CovertLab/borealis) for the code and documentation
on running FireWorks workflows in Google Cloud including **Team Setup** and
**Developer Setup** steps.


## One time setup

1. [Install the Google Cloud SDK and log in](https://github.com/CovertLab/borealis/blob/master/docs/install-tools.md).
   The SDK includes the `gcloud` and `gsutil` command line programs.
   `gcloud` commands configure and operate projects, access controls, and servers.
   `gsutil` commands manipulate Google Cloud Storage (GCS) buckets, transfer files,
   run rysnc, etc.
   * See the note in that webpage about working around a compatibility bug by
     configuring `CLOUDSDK_PYTHON` to Python 3.7.
   * `gcloud` will occasionally prompt to install updates.
     You can install updates proactively via:  
     `gcloud components update`
   * The [Cloud Console](https://console.cloud.google.com/home/dashboard) is a web
     interface to manage your cloud services, view logs, etc.
     Open the ☰ (hamburger) menu in the top left corner to access all the
     sub-console pages.
1. Create your own Google Cloud Storage bucket for output files.
   Using a bucket per developer enables usage tracking, cleanup, and ACLs.
   *  You can create it via the [Google Cloud Storage — Browser
      webpage](https://console.cloud.google.com/storage/browser) or via the
      `gsutil` command line program.
   *  Pick a name like "sisyphus-crick".  
      It has to be globally unique.
      BEWARE: The name is **publicly visible** so don't
      include login IDs, email addresses, project names, project numbers, or
      personally identifiable information (PII).
   *  Pick the same Region used with Compute Engine (run `gcloud info` for info).
   *  Pick Standard storage class, and default access control.
   *  Store the name in the environment variable `WORKFLOW_STORAGE_ROOT` in
      your shell profile and update your current shell, e.g.:  
      `export WORKFLOW_STORAGE_ROOT="sisyphus-crick"`
   *  If you don't have access permissions to create a storage bucket, ask a gcloud project
      administrator. As a fallback, you can use a subdirectory of the `sisyphus`
      bucket, e.g. `sisyphus/data/crick`.


## Setup

* Before doing the **Run it** steps in the next section, set your current working directory to
  your wcEcoli git clone, e.g.:  
  `cd ~/dev/wcEcoli`
* Running Python in the next section needs the wcEcoli directory on the Python path:  
  `export PYTHONPATH=$PWD`  
  A shell alias is useful for this:  
  `alias ppath='export PYTHONPATH=$PWD'`
* Create a `my_launchpad.yaml` file like this to use with our MongoDB LaunchPad in
  Google Compute Engine:

  ```
  host: localhost
  name: YOUR-USERNAME-HERE
  username: null
  password: null
  port: 27017
  strm_lvl: INFO
  ```

  This has `host: localhost` for use with the `mongo-ssh.sh` script, below, which opens
  an SSH tunnel from `localhost:27017` to port `27017` on the MongoDB server.

  (A LaunchPad database can support multiple workflows although the worker nodes will
  fetch Firetasks from all the READY tasks in that database. If you want another
  LaunchPad for an independent set of worker nodes, create another launchpad yaml file
  and switch between launchpads via command line options to the `lpad`, `wcm.py`, and
  `gce` commands.)


## Run it

Also see [Borealis: How to run a workflow](https://github.com/CovertLab/borealis#how-to-run-a-workflow).

**Note:** The script `cloud/build.sh` runs both step 1 and step 2.

1. Build your `"$USER-wcm-runtime"` Docker Image:

   ```shell script
   cloud/build-runtime.sh
   ```

   This builds a Docker Image containing Python, binary libraries, and pip packages,
   building the equivalent of
   [Creating the "pyenv" runtime environment](create-pyenv.md) for worker nodes.  

   It does not contain your Python source code in the project.  
   This step takes about 20 minutes, but the heavy lifting happens in a Google server.

   Rerun `cloud/build-runtime.sh` whenever the Python executable or libraries change (i.e. when
   you'd update your local pyenv virtualenv) and you want to run it in Google Cloud.
   To use a new `wcm-runtime`, you'll need to build a new `wcm-code` on top of it.

   This builds a fresh Docker Image as defined by the files
   `cloud/docker/runtime/Dockerfile` and `requirements.txt` in your current git
   working directory.

1. Build your `"$USER-wcm-code"` Docker Image to deploy your Whole Cell Model code:

   ```shell script
   cloud/build-wcm.sh
   ```

   This builds a Docker Image starting from `"$USER-wcm-runtime"`, adding the
   Python source code in your current git working directory and compiling its Cython code.  
   This step only takes a few minutes.

   (Your code does not need to be checked into git but consider testing it with `pytest`,
   `mypy`, and some manual runscripts first. You can quickly run a partial
   sim generation via  
   `python runscripts/manual/runSim.py --length-sec 120 --no-jit`)

   Rerun `cloud/build-wcm.sh` whenever the `"$USER-wcm-runtime"` Image changes or your
   Python source code changes that you want to run in Google Cloud.  
   Cloud Firetasks that start running after you build an Image will "pull" it and use it,
   even Firetasks that were already queued to run. You can do this to update your code
   then use `lpad` commands to restart failed tasks and resume paused tasks.  
   The full Docker Image path is `gcr.io/$PROJECT/$USER-wcm-code`, where
   `$PROJECT` is the Google Cloud project name and `$USER` is your local username.
   You can also "pull" this Image and run it locally.

1. Open an ssh tunnel to the project's MongoDB "LaunchPad" server in Google Compute Engine:

   ```shell script
   runscripts/cloud/mongo-ssh.sh
   ```

   You only need this connection while you're uploading Fireworks workflows,
   examining them with commands like `lpad get_fws`, viewing workflow state via
   `lpad webgui`, or otherwise accessing the LaunchPad. It's fine to leave this
   connection open for days. It's also fine to close this connection while a
   workflow continues to run.

   **Tip:** You can run `runscripts/cloud/mongo-ssh.sh bg` to open the connection
   and leave it open in a background process without tying up a shell tab.

1. Initialize your MongoDB database before first use:

   ```shell script
   lpad reset
   ```

   The simplest way to use Fireworks is to `reset` the LaunchPad like this before
   uploading a workflow (instructions below). But there are more patterns of
   Fireworks use,
   e.g. you can let workflows accumulate in the database, rerun them later, or
   delete or archive specific workflows.

   **Tip:** If you `reset` often enough, you might like this way to breeze
   past its confirmation step:

   ```shell script
   echo y | lpad reset
   ```

1. Upload and start the wcEcoli workflow from another terminal tab:

   ```shell script
   python runscripts/cloud/wcm.py
   ```

   * Use the `-h` option to get help on all the command line options, e.g. number
   of `--init-sims` initial lineage seeds to run, number of simulation `--generations`
   to run per seed, and a workflow `--description`.
   * The `--workers` option sets the number of worker nodes to launch to run this
   workflow, otherwise `wcm.py` will use a heuristic. Use `--workers 0` (or `-w0`)
   if you don't want to launch worker nodes yet.
   * Use the `--dump` option to dump the workflow definition as a YAML file for
   review _instead of_ uploading it to the LaunchPad.
   * After `wcm.py` uploads the workflow, you can close the `mongo-ssh.sh` ssh tunnel
   or leave it open for interactive commands like `lpad get_fws`.
   * Use the `gce` command to launch (more) workers when you want.  
     **Tip:** Launching workers via the `gce` command is especially
     useful after a task failure if the workers timed out and shut down by the time
     you upload a fixed `"$USER-wcm-runtime"` Docker Image.


## Watch it run

### Setup

* In your web browser, bookmark the [Google Cloud Platform
console](https://console.cloud.google.com/home/dashboard)
and log in.
* In the `☰` (hamburger) menu in the top left corner of the Google Cloud Platform console webpage,
pin "Compute Engine" (GCE), "Logging" (Stackdriver), and "Storage" (GCS) to the top section of this menu
for quick access.
* You can rearrange and hide the "cards" on the console's home page.

**Optional:** Set up a chart to display CPU utilization.

  1. Move the "Compute Engine" card to top center position.
  1. In this card's `⋮` menu, click `Add Chart`.
  1. Title it `CPU usage`.
  1. Add Metric type `instance/cpu/utilization` with
     filter `resource.type = "gce_instance"`.
  1. Click `Save`.

To make a chart for just Fireworkers, user a filter like
`metric.labels.instance_name = starts_with("fireworker")`.


### Monitor

* Use your web browser bookmark to open the [Google Cloud Platform
console](https://console.cloud.google.com/home/dashboard)
home page. Use the `☰` menu or more bookmarks to navigate to the other pages.

* Open the [Compute Engine — VM instances](https://console.cloud.google.com/compute/instances)
page to see the list of running Compute Engine VM instances.

   * The `fireworker-USER-0`, `fireworker-USER-1`, ... VM instances are the worker
   nodes to run the workflow firetasks, where `USER`
   will be the `name:` value in your `my_launchpad.yaml` file, typically set to
   your shell `$USER` username.
   Use `wcm.py` or `gce.py` to launch worker nodes.
   * "mongo-prime" is running MongoDB for the FireWorks LaunchPad.

* Open the [Google Cloud Logs Explorer](https://console.cloud.google.com/logs/)
page to view the logs from the project's GCE VM instances.

   * The "Stream logs" button at the top starts streaming the logs (akin to `tail -f` following).
   * See [Debugging](#Debugging) for tips on filtering the logs.
   * Each firetask writes a log file to the `logs/` output subdirectory. See [Debugging](#Debugging).

* Open the [Storage — Browser](https://console.cloud.google.com/storage/browser)
page to browse the files created by the workflow.

   * `wcm.py` directs the workflow's output files to the storage bucket you picked, above,
   and a subdirectory named `WCM/$TIMESTAMP/` or `WCM/$TIMESTAMP__description/`.
   * Its `logs/` subdirectory contains a log file for each Firetask run, whether
   successful or not, but each log file doesn't get stored until the Firetask ends.


## Download the outputs

### Ways to download the workflow outputs

* **Simplest for individual files:** Open the [Google Cloud
Storage — Browser](https://console.cloud.google.com/storage/browser),
browse into your Google Cloud Storage bucket ("sisyphus-crick" or whatever)
find your workflow files, and click on individual files to download them.
* **Most convenient:** Use [gcsfuse](https://github.com/GoogleCloudPlatform/gcsfuse) to mount your
storage bucket to your local file system (or use an option like `--only-dir WCM`
to mount just a subdirectory of it) and access
the files like local files.
  * Example:  
    `cd ~/dev/gcs/ && mkdir sisyphus-crick && gcsfuse sisyphus-crick sisyphus-crick`
  * Note the **gcsfuse setup** instructions below.
  * gcsfuse can't completely hide the fact that Google Cloud Storage (GCS) is not a
  regular file system.
  E.g. GCS doesn't have directories, just file paths that may contain slashes and may end with
  a slash, so it cannot atomically rename a "directory".
  GCS reads and writes whole files.
  See the notes in
  [semantics.md](https://github.com/GoogleCloudPlatform/gcsfuse/blob/master/docs/semantics.md)
  for details.
* **Fastest bulk transfers:** Use the
[gsutil command line tool](https://cloud.google.com/storage/docs/gsutil). (It's part of the
Cloud SDK tools you installed along with `gcloud`.)
  * Use `gsutil -m` to transfer multiple files in parallel.
    Example: `gsutil -m cp gs://my-bucket/my-subdir download_dir`
  * It can also do [gsutil -m rsync](https://cloud.google.com/storage/docs/gsutil/commands/rsync).
* We could write a visual browser for images in our Cloud Storage bucket.
* If you have lots of microscopy images, try running [ffmpeg](https://ffmpeg.org/) to
compress them as as video using a lossless or "visually lossless" codec such as
[LosslessH.264](https://trac.ffmpeg.org/wiki/Encode/H.264#LosslessH.264),
[FFV1](https://trac.ffmpeg.org/wiki/Encode/FFV1), or
[others](http://compression.ru/video/codec_comparison/lossless_codecs_2007_en.html)
to save on storage and transfer anywhere.


### gcsfuse setup

See their [docs to install it](https://github.com/GoogleCloudPlatform/gcsfuse/blob/master/docs/installing.md).

Authenticating with **gcsfuse** requires a credentials key file for a Google Cloud
_service account._ We use the `fireworker` service account since it can not only access
Google Cloud Storage, it has the Logs permissions needed to run a fireworker, in case
you want to do that locally.

Run these shell commands to get a credentials key file `fireworker.json`:

```shell script
export GOOGLE_APPLICATION_CREDENTIALS=$HOME/bin/fireworker.json
gcloud iam service-accounts keys create $GOOGLE_APPLICATION_CREDENTIALS --iam-account=fireworker@allen-discovery-center-mcovert.iam.gserviceaccount.com
chmod 400 $GOOGLE_APPLICATION_CREDENTIALS  # owner-only since it contains a private key
# Now add that "export" line to your shell .profile or .bash_profile file.
```

(_Alternatively_ you can get the credentials key file from the Google Cloud Console's
IAM-Service Accounts section: Click on the "fireworker" service account,
click `ADD KEY` > `Create new key` > `JSON` > `CREATE`, and put the downloaded file at
`$HOME/bin/fireworker.json`. As above: `chmod` it, export/set
`$GOOGLE_APPLICATION_CREDENTIALS`, and add that "export" line to your shell profile.)


## Debugging

Use the [Logs Explorer](https://console.cloud.google.com/logs/) to watch workflows
run with error messages and other information from the servers involved.

(Set a Logs query like `resource.labels.instance_id:"fireworker" severity>=INFO` to
see just the big picture of Fireworkers and Firetasks starting and stopping.
The `:` means "starts with."  
Expand the log level to `severity>=DEBUG` to include task stdout and stderr lines.  
If coworkers are also running workflows, you can filter to just your own
Fireworkers using a query like
`resource.labels.instance_id:"fireworker-USER-"`,
inserting your $USER name.)

See the [Download the outputs](#Download-the-outputs) for ways to download the
`logs/*.log` files of the firetask  runs. If you rerun firetasks (e.g.
`lpad rerun_fws -i 10,12`), each run will write another log file.

After mounting the GCS storage bucket via **gcsfuse** or downloading its `logs/` files,
you can quickly check on completed firetask runs:

```shell script
grep TASK: logs/*           # see SUCCESSFUL and FAILED runs
grep -h TASK: logs/*        # suppress the filenames for easier scanning through them
grep '[A-Z]* TASK:' logs/*  # make the "SUCCESSFUL" and "FAILED" part stand out more
grep FAILED logs/*          # see only the FAILED runs
```

Use FireWorks `lpad` commands to view the state of your Firetasks. See
[Borealis: How to run a workflow](https://github.com/CovertLab/borealis#how-to-run-a-workflow).

To change your code and run again:

1. Rerun `cloud/build-wcm.sh` to push a new Docker Image.
1. Mark the relevant Firetasks to run again, for example:
   ```shell script
   lpad rerun_fws -i <FW_IDS>
   ```
   or
   ```shell script
   lpad rerun_fws -s FIZZLED
   ```
1. Check the [Compute Engine VM instances](https://console.cloud.google.com/compute/instances)
or run the command
   ```shell script
   gcloud compute instances list
   ```
to see if the worker nodes are still running. If not, use the Borealis `gce` command
to launch some workers.


## Notes

When using a Docker image name with no :tag, Docker applies the default
tag "`:latest`".
**[It does not mean the "latest version"](https://vsupalov.com/docker-latest-tag/)!**
It's just a default tag name. Pulling with the default tag or with `:latest` (same
thing) will pull the latest image that has the tag `:latest`, as with any tag.
It won't pull any image that has a different tag.


### Multiple `wcm-code` Docker images

If you want to build more than one `wcm-code` Docker image to run different code in
different workflows, pass an optional "workflow user ID" to `cloud/build-wcm.sh`, e.g.:

   ```shell script
   cloud/build-wcm.sh nightly
   ```

and pass that workflow user ID to `wcm.py` via the environment variable `$WF_ID`:

   ```shell script
   WF_ID=nightly python runscripts/cloud/wcm.py
   ```
