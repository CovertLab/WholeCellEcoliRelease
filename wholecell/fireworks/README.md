# One-time setup

## MongoDB

1. Create an account at mongolab.com

2. Sign in and get to the home screen

3. Next to "MongoDB Deployments" you'll see three buttons. Click the one that says "Create new".

4. For "Cloud provider", select "amazon web services".  For "Location": "Amazon's US East (Virginia) Region (us-east-1)".

5. Under "Plan", choose "Single-node" and select "Sandbox" (...it's free).

6. For "Database name" write "wc_ecoli".

7. Click "Create new MongoDB deployment".

8. Back on the home screen, click on "wc_ecoli".

9. Click on the "Users" tab and then select "Add a database user".

10. Choose a username and password.  Note that in fireworks, the password is stored in plaintext.

11. Note that the information shown at the top ("To connect using the shell") contains the database hostname and database port (the part after the colon is the port).

12. Make a SECOND MongoDB Deployment, by repeating steps 3-10. Everything is the same, except for step 6 - give the new deployment the name "wc_ecoli_2"

## Config files for Fireworks

* Run `python wholecell/fireworks/initialize.py`

* This will ask for two sets of data - input numbers for the first deployment you made first, and for the second deployment you made second

* Run `lpad -l my_launchpad.yaml reset` and choose `Y`es if prompted

* Run `lpad -l my_launchpad_2.yaml reset` and again choose `Y`es if prompted