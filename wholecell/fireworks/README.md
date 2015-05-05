# One-time setup

## MongoDB

* Create an account at mongolab.com

* Sign in and get to the home screen

* Next to "MongoDB Deployments" you'll see three buttons. Click the one that says "Create new".

* For "Cloud provider", select "amazon web services".  For "Location": "Amazon's US East (Virginia) Region (us-east-1)".

* Under "Plan", choose "Single-node" and select "Sandbox" (...it's free).

* For "Database name" write "wc_ecoli".

* Click "Create new MongoDB deployment".

* Back on the home screen, click on "wc_ecoli".

* Click on the "Users" tab and then select "Add a database user".

* Choose a username and password.  Note that in fireworks, the password is stored in plaintext.

* Note that the information shown at the top ("To connect using the shell") contains the database hostname and database port (the part after the colon is the port).

## Config files for Fireworks

* Run `python wholecell/fireworks/initialize.py`

* Run `lpad reset` and choose `Y`es if prompted