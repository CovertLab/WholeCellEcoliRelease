# One-time setup

To use FireWorks, your computer (to run the `lpad` command) and
the Firetask worker nodes need to access a MongoDB
database holding the workflow information.


## Creating a MongoDB database

The team has a MongoDB server in the Allen Center project on Google Cloud.
It's useful for running FireWorks in Google Cloud and on our local computers.
(See [How to run the Whole Cell Model on the Google Cloud
Platform](../../docs/google-cloud.md) for details.)
You can access the service securely via an ssh tunnel with your Google login.
Creating a database there just takes the steps in the next section.

That MongoDB service is not accessible to the open Intenet,
so to run workflows on Sherlock you can create an account and database on
[MongoDB Atlas](https://www.mongodb.com/cloud/atlas).
See the instructions in the
[MAT project forum](https://matsci.org/t/heres-how-to-connect-to-atlas-mongodb/4816)
to set it up including using the `lpad init -u` command to create a
`my_launchpad.yaml` that configures FireWorks to connect to it.

**CAUTION:** Your Atlas database will be reachable to the entire Internet.
So use a password vault to generate and save long, random passwords for the
account and the DB user. As always, don't reuse them anywhere else.


## Config YAML files for Fireworks

Fireworks needs a launchpad config file, `my_launchpad.yaml`.
When running with a job queue via `qlaunch` (such as with SLURM on Sherlock, not
with worker nodes on Google Compute Engine),
it also needs a queue-adapter file, `my_qadapter.yaml`.

1. Create the `my_launchpad.yaml` file.

   * **Atlas MongoDB cluster:** See the
     [MAT project forum](https://matsci.org/t/heres-how-to-connect-to-atlas-mongodb/4816)
     on how to use the `lpad init -u` command to create `my_launchpad.yaml`.

   * **Google Cloud MongoDB database:** Run the script
     `runscripts/cloud/mongo-ssh.sh` to open an ssh tunnel whenever you want to
     access the server. (You can leave it running in one terminal tab then `^C`
     it when you're done, or run `runscripts/cloud/mongo-ssh.sh bg` to run it in a
     background process and kill it when you're done.)

     Construct `my_launchpad.yaml` like this:

     ```
     host: localhost
     port: 27017
     name: <your-username-or-another-unique-name-for-your-database>
     username: null
     password: null
     ```

     These values ask FireWorks to connect to `localhost:27017`, which the
     `mongo-ssh.sh` tunnel will forward to port `27017` on the Google Compute
     Engine MongoDB server.

     `name` will be your database name.

     We're not using user
     authentication on this server since the ssh tunnel already requires Google login.

1. Run

   ```
   lpad reset
   ```

   This will connect to the MongoDB server then prompt you to type `Y`es to
   confirm creating and initializing the database. If it fails to connect, check
   that the `mongo-ssh.sh` tunnel is still running.

   If you want additional launchpad databases, create a launchpad yaml file
   for each one.
   To use a launchpad config filename besides the default
   `my_launchpad.yaml`, pass it as a command option like
   `lpad -l gce_launchpad.yaml webgui`.

1. If you need a queue-adapter to run FireWorks workers on Sherlock, write
   `my_qadapter.yaml` from the template
   `wholecell/fireworks/templates/my_qadapter.yaml`.
   Or you can run the interactive script `python -m wholecell.fireworks.initialize`,
   but beware that it overwrites `my_launchpad.yaml` and `my_qadapter.yaml`.

The launchpad database keeps the status of current and past workflows.
You can rerun past workflows without re-uploading them, archive or delete
workflows, and mix multiple workflows in the same database, but the simplest
approach is to `lpad reset` it and upload a new workflow each time.
