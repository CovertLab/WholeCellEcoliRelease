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

In your cloned `wcEcoli` directory, to compile cython plugins and any C functions with python bindings, run:

```bash
make compile
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



Running a simulation on a second launchpad
------

If you made a second launchpad while setting up fireworks (which you did if you followed the README.md in the fireworks folder),
then it's possible to queue and run from two launchpads at once. For example this would allow running one long simulation over
the course of several days while also running others to work on code and debug in the meantime.

To check if you set up a second launchpad, look in your wcEcoli folder. There should be a my_launchpad.yaml and a my_qadapter.yaml file.
If there are also my_launchpad_2.yaml and my_qadapter_2.yaml files, then you have configured a second launchpad.

To queue a task into your second launchpad, use the LAUNCHPAD_FILE flag when you call fw_queue.py:

```bash
DESC="Example run on a second launchpad." LAUNCHPAD_FILE=my_launchpad_2.yaml python runscripts/fw_queue.py
```

To run any lpad commands on the second launchpad (such as lpad get_fws, or lpad reset), use the -l flag and add the name of the second launchpad yaml file:

```bash
lpad -l my_launchpad_2.yaml get_fws
```

If you run commands without specifying, it will default to using my_launchpad.yaml, which is your primary launchpad.

To rlaunch from a specific launchpad, use the same -l flag:

```bash
rlaunch -l my_launchpad_2.yaml rapidfire
```

To qlaunch from a specific launchpad, again use -l and the name of the yaml file corresponding to that launchpad. For this command to work,:

```bash
qlaunch -l my_launchpad_2.yaml -r rapidfire --nlaunches infinite --sleep 5
```
NOTE: while this will give you a separate fireworks launchpad, things launched from that launchpad will still use the code in your
directory, even if you've modified it since you prepared the fireworks queue. In other words, if you are running a long set of 
simulations and debugging while you do so, then later generations of the long simulation will start running on your modified code.

One solution to this problem is to clone a separate repo fresh from a branch and launch your alternate workflow from that location,
so that you will not be updating the code on which your simulation will run.