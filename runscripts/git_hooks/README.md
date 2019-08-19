# Git hooks
This directory contains useful git hooks for developing in the repo on Linux and Mac systems.  Copy the scripts in this directory to your git hooks directory (likely `wcEcoli/.git/hooks/`).  Ensure that the hooks have execution permissions.

From the `wcEcoli/` directory, run the following command to copy the scripts:
```
cp runscripts/git_hooks/*[^.md] .git/hooks/
```

## Hooks

* `post-checkout`: Runs when performing `git checkout` to a different branch
* `post-merge`: Runs when merging commits

## Scripts
* `wcm-env-check.sh`: Called by `post-checkout` and `post-merge` to ensure an up to date environment for development.  Does the following:
    * Remove any `*.pyc` files to prevent python executing these even if the corresponding `*.py` file has been removed.
    * `make compile`: ensure compiled code is up to date if changed
    * Show difference between packages in `requirements.txt` and those that are pip installed. Try pip installing missing or out of date packages if any are displayed.
