# One-time setup

## Creating a MongoDB database on mlab

**NOTE:** We've been using databases on mlab.com for running FireWorks on Sherlock.
We also have MongoDB running in the Allen Center project on Google Compute
Engine. With that server, just pick a unique database `name` such as your username,
skip over these mlab steps, and create the launchpad config file (down below).


1. Create an account at [mlab.com](https://mlab.com/).

2. Sign in.

3. Next to **MongoDB Deployments** page, click the `+ Create new` button.
   Do _not_ click the `+ Create new` button next to **Environments**.

4. Pick a "Cloud provider". We've been using `amazon web services`.

5. Pick the "Plan Type" `SANDBOX [FREE]`.
   It shows which regions are available for this Plan Type.

   Then click the blue `CONTINUE` button way down at the bottom of the page.

7. Pick a "Region". We've been using Amazon's `US East (Virginia) (us-east-1)`.

   Then click the blue `CONTINUE` button again at the bottom of the page.

8. Pick a "Database name" such as `wc_ecoli`.

   Then click the blue `CONTINUE` button again at the bottom of the page.

   It should confirm the Plan Type, Cloud Provider, Region, Plan size 0.5 GB
   (which is plenty), MongoDB version, database name, and price **FREE**.

9. Click `SUBMIT ORDER`.

   Returning you to the "MongoDB Deployments" page, it should show the new
database in the list and it should (soon) show "Ok: This database is up and running".

10. Click on the database, then click the "Users" tab, then click "Add database user".

11. Choose a username (e.g. `fireworks`) and a new, long, random password.

    **CAUTION:** This database will be reachable to the entire Internet and the
    fireworks configuration files you'll create in a moment will store the
    database password in plain text. So use a password vault to generate a new,
    long, random password and (as always) don't reuse it anywhere else.

12. Copy the information shown at the top of the page ("To connect ...") with
the database hostname and database port (the number after the ":" is the port)
to your password vault, along with the database name, database username, and
database user password.

13. Optional: If you want to be able to run two workflows at the same time, make a
second MongoDB Deployment by repeating the new-deplohment steps above with a
different database name such as `wc_ecoli_2`.


## Launchpad config YAML files for Fireworks

Fireworks needs a launchpad config file, by default named `my_launchpad.yaml`.
When running with a job queue via `qlaunch` (such as with SLURM on Sherlock),
it also needs a queueadapter file, by default named `my_qadapter.yaml`. The
queueadapter is not currently needed when running on Google Compute Engine (GCE).

1. Either run the interactive script `python -m wholecell.fireworks.initialize`
to construct `my_launchpad.yaml` and `my_qadapter.yaml` files. It will
overwrite existing files. _Or_ you can manually crib from the template files
`wholecell/fireworks/templates/*.yaml`.

   Running `initialize.py`, enter values when prompted. The default values are
   suitable for using MongoDB on `localhost` or connecting to it on
   Google Compute Engine via `mongo-ssh.sh` ssh port forwarding tunnel.
   The database name defaults to `$USER` for a user-specific database.
   Or enter the mlab details if you're using mlab.com.

   To run on GCE, you just need a LaunchPad YAML file (the default name is
   `my_launchpad.yaml`) containing:

       host: localhost
       name: <your-username-or-another-unique-name>
       username: null
       password: null
       port: 27017

2. Run `lpad reset`. This will connect to the database then
prompt you to choose `Y`es to confirm resetting the workflow state.

   If you want multiple MongoDB launchpad databases, e.g. for independent
   FireWorks workflows, create a launchpad yaml file for each one. On GCE
   you just need to pick a unique `name:` for each one. You can crib from your
   first `my_launchpad.yaml` file or re-run `initialize.py`, but beware that
   `initialize.py` will overwrite `my_launchpad.yaml` and `my_qadapter.yaml`.

   To use a launchpad config file besides the default filename
   `my_launchpad.yaml`, pass it as an option to commands, e.g.
   `lpad -l gce_launchpad.yaml reset`.

   FYI, a single launchpad database can rerun past workflows and support multiple
   workflows, but the simplest approach is to `reset` it and upload a new workflow
   each time.
