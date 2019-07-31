# One-time setup

## MongoDB

1. Create an account at [mlab.com](https://mlab.com/).

2. Sign in.

3. On the "MongoDB Deployments" page, click the `+ Create new` button.

4. For "Cloud provider", select `amazon web services`.

5. For "Plan Type", choose `SANDBOX [FREE]`.

6. Click the blue `CONTINUE` button way down at the bottom of the page.

7. Pick a location. We usually pick Amazon's `US East (Virginia) (us-east-1)`
then click the blue `CONTINUE` button again at the bottom of the page.

8. For "Database name" enter `wc_ecoli`
then click the blue `CONTINUE` button again at the bottom of the page.

   It should confirm the Plan Type, Cloud Provider, Region, Plan size 0.5 GB
   (which is plenty), MongoDB version, database name, and price FREE.

9. Click `SUBMIT ORDER`.

   Returning you to the "MongoDB Deployments" page, it should show the new
database in the list and it should (soon) show "Ok: This database is up and running".

10. Click on the database, then click the "Users" tab, then click "Add database user".

11. Choose a username (e.g. `fireworks`) and a new, random password.

    **CAUTION:** The fireworks configuration files you'll create in a moment will store
    the database password in plain text. So use a password vault to generate a new,
    random password and (as always) don't reuse it anywhere else.

12. Copy the information shown at the top ("To connect to mongo shell" and
"To connect using a driver via the standard MongoDB URI") that contain the database
hostname and database port. (The number after the ":" is the port.)

13. Optional: If you want to be able to run two workflows at the same time, make a
second MongoDB Deployment by repeating steps 3-12 with a different database name
`wc_ecoli_2` in step 8.

## Config files for Fireworks

* Run `python wholecell/fireworks/initialize.py`

* Run `lpad -l my_launchpad.yaml reset`. This will connect to the database then
prompt you to choose `Y`es to confirm resetting the workflow state.

* Optional: Run `lpad -l my_launchpad_2.yaml reset` and again choose `Y`es if prompted
