Usage
======

To run simulations on Sherlock, do the following:

One time setup
--------------
In your cloned `wcEcoli` directory, to set up the proper python environment, run:

```bash
pyenv local wcEcoli
```

Complete the one-time setup for fireworks as described in [wholecell/fireworks/README.md](wholecell/fireworks/README.md)

You will also need to add the following to your `$HOME/.bash_profile` (using the appropriate path):

```bash
export PYTHONPATH="/path/to/wcEcoli:$PYTHONPATH"
```

To import the necessary shared libraries, you will need to execute the following *each time* after you log in to Sherlock (alternatively, you can add it as a line to your `$HOME/.bash_profile`):

```bash
module load wcEcoli
```

Single simulation
------------------

To queue a simulation in fireworks:

```bash
DESC="Example run of a single simulation." python runscripts/fw_queue.py
```

Note that the text provided to the `DESC` variable should be changed to something more descriptive.

Multiple simulations
--------------------

To queue multiple simulations (in this case 4 simulations) in fireworks:

```bash
DESC="Example run of multiple simulations." N_INIT_SIMS=4 python runscripts/fw_queue.py
```

Multiple generations
--------------------

To queue multiple generations (in this case 4 generations) from a single mother cell:

```bash
DESC="Example run of multiple generations." N_GENS=4 python runscripts/fw_queue.py
```

To queue multiple generations (in this case 3 generations) from multiple mother cells (in this case 2 mother cells:

```bash
DESC="Example run of multiple generations from multiple mother cells." N_GENS=3 N_INIT_SIMS=2 python runscripts/fw_queue.py
```


Using an interactive node to run simulations
--------------------------------------------

To run simulations on an interactive session (after having queued them), run:

```bash
rlaunch rapidfire
```

You probably only want to do this if you're running/debugging a single simulation.

Don't do this on a login node.

Using the scheduler (SLURM) to run simulations
-----------------------------------------------

To run simulations using the cluster (you'll probably want to do this if you're running more than one simulation and/or more than one generation), run:

```bash
qlaunch -r rapidfire --nlaunches infinite --sleep 5
```

This command will run forever until you `Ctrl-C` to kill it once you see that all the output and analysis files have been generated.

`qlaunch` is relatively lightweight, so you can probably get away with running it on a login node.


Output
------

The output is stored as a time-stamped sub-directory of the `out` directory, for example `out/20140825.095758.954820584`.

Within this directory, there is a `metadata` sub-directory which stores the git revision information as well as the description provided by the `DESC` variable, a `kb` sub-directory which stores kb objects (after the simulations and analysis are done the objects are compressed using bzip2), and sub-directories (maybe only a single sub-directory) containing different variants (e.g., gene knockouts or other perturbations).

Within variant sub-directories, there are `N_INIT_SIMS` (which defaults to 1) numbered sub-directories such as `000000` corresponding to "family trees".

Within each "family tree" sub-directory are generation sub-directories such as `generation_000000`.

Within each generation sub-directory are numbered individual simulation directories that contain `simOut` (for simulation data) and `plotOut` (for plots) sub-directories.

A family tree for 3 generations showing the relationship between numbered individual simulations is shown here:

```
gen 0:                 0
gen 1:         0               1
gen 2:     0       1       2       3
```