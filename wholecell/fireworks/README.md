# One-time setup

## Make a MongoDB database

Running within FireWorks requires a MongoDB database. You can install a MongoDB
server, use a shared server, or instantiate a database on mlab.com.

1. Create an account at mlab.com

2. Sign in and get to the home screen

3. Next to "Create New Deployment".

4. For "Cloud provider", select "amazon web services".  For "Location": "Amazon's US East (Virginia) Region (us-east-1)".

5. Under "Plan", choose "Sandbox" (...it's free).

6. For "Database name" write "wcecoli".

7. Click "Create new".

8. Back on the home screen, click on "wcecoli".

9. Click on the "Users" tab and then select "Add a database user".

10. Choose a username and password.  Note that in FireWorks, the password is stored in plaintext.

11. Note that the information shown at the top ("To connect to mongo shell") contains the database hostname and database port (the part after the colon is the port).

12. Optional: Make a SECOND MongoDB Deployment, by repeating steps 3-10. Everything is the same, except for step 6 - give the new deployment the name "wc_ecoli_2".  This is useful if running from a second repo and not necessary with only one clone.

## Make the Configuration files for FireWorks

* Run `python wholecell/fireworks/initialize.py`

* Run `lpad -l my_launchpad.yaml reset` and choose `Y`es when prompted
