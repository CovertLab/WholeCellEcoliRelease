Whole Cell Model - *Escherichia coli*
===================================

This repository contains work to date on the [Covert Lab's](https://www.covert.stanford.edu/) Whole Cell Model for [*Escherichia coli*](https://en.wikipedia.org/wiki/Escherichia_coli), as well as some effort to create a framework for building whole cell models in general. 

Until now most of the work has been on Stanford's local [Sherlock](https://www.sherlock.stanford.edu/) cluster, though recent efforts have enabled the simulation and its associated analyses to be run manually on a local machine.

Usage
======

To run simulations on Sherlock, do the following:

One time setup
--------------
In your cloned `wcEcoli` directory, to set up the proper python environment, run:

```bash
pyenv local wcEcoli2
```

Complete the one-time setup for fireworks as described in [wholecell/fireworks/README.md](wholecell/fireworks/README.md)

To access or install the necessary shared libraries...
  * on Sherlock, execute the following *each time* you log in to Sherlock (*tip:* add it to your `$HOME/.bash_profile`):

	```bash
	module load wcEcoli/sherlock2
	```

  * on another computer, see [How to set up the runtime environment for the model](https://github.com/CovertLab/wcEcoli/wiki/How-to-set-up-the-runtime-environment-for-the-model) in the wiki.

Set up the $PYTHONPATH like this (*tip:* add this to your `$HOME/.bash_profile`):

```bash
export PYTHONPATH="/path/to/wcEcoli:$PYTHONPATH"
```

In your cloned `wcEcoli` directory, to compile cython plugins and any C functions with python bindings, run:

```bash
make compile
```

You need to make sure the output from the model goes to the SCRATCH filesystem (which is larger) rather than SHERLOCK HOME. You'll need to make a symbolic link between the output directory of your `wcEcoli/` directory and a directory in SCRATCH. Within your wcEcoli diretory, there needs to be a folder named `out/`. The model puts its output in this folder, so we can basically tell the computer to send anything placed in this folder to SCRATCH instead. You can call the folder on SCRATCH whatever you want, but here is one option:

```bash
mkdir $SCRATCH/wcEcoli_out
```

Now, make sure an out directory in the cloned wcEcoli directory doesn't exist already and create a new symbolic link from out to the directory on SCRATCH that was just created above:

```bash
ln -s $SCRATCH/wcEcoli_out out
```

Similarly, we would like to create a symbolic link to a shared sim data cache directory on PI_SCRATCH that should contain a copy of the newest sim data object (it should be updated daily with any new changes to the codebase):

```bash
ln -s $PI_SCRATCH/wc_ecoli/cached cached
```

Running a simulation
-------------------- 

There are two ways to run a simulation:

   1. Queue up a Fireworks workflow, then run it.
   2. Use the manual runscripts.

The command line program `fw_queue.py` queues up a Fireworks workflow including parameter fitting, the simulation itself, and lots of analysis plots.

The source file `fw_queue.py` begins with documentation on its many options. Below are a few usage examples.

But first, note that you can reset the Fireworks queue (if needed) via:

```bash
lpad reset
```

### Single simulation

To queue up a single simulation in Fireworks:

```bash
DESC="Example run of a single simulation." python runscripts/fireworks/fw_queue.py
```

The `DESC` text should be more descriptive so you can readily distinguish your runs.

### Multiple simulations

To queue multiple simulations (in this case 4 simulations) in fireworks:

```bash
DESC="Example run of multiple simulations." N_INIT_SIMS=4 python runscripts/fireworks/fw_queue.py
```

### Multiple generations

To queue multiple generations (in this case 4 generations) from a single mother cell:

```bash
DESC="Example run of multiple generations." N_GENS=4 python runscripts/fireworks/fw_queue.py
```

To queue multiple generations (in this case 3 generations) from multiple mother cells (in this case 2 mother cells:

```bash
DESC="Example run of multiple generations from multiple mother cells." N_GENS=3 N_INIT_SIMS=2 python runscripts/fireworks/fw_queue.py
```

### Shifting nutrient conditions

To queue a simulation that switches between environments, use the "nutrient_time_series" variant, and give the range of indices (in this case from 1 to 1) specifying conditions defined in wcEcoli/reconstruction/ecoli/flat/condition/timeseries:

```bash
DESC="Example run of nutrient shifts." VARIANT="nutrient_time_series" FIRST_VARIANT_INDEX=1 LAST_VARIANT_INDEX=1 python runscripts/fireworks/fw_queue.py
```

### Using the cached sim data object

To use the cached sim data object, use the CACHED_SIM_DATA environment variable:

```bash
DESC="Example run with cached sim data." CACHED_SIM_DATA=1 python runscripts/fireworks/fw_queue.py
```

### Using an interactive node to run a Fireworks workflow

To run simulations on an interactive session (after having queued them), run:

```bash
rlaunch rapidfire
```

You probably only want to do this if you're running/debugging a single simulation.

Don't do this on a Sherlock login node.

### Using the scheduler (SLURM) to run a Fireworks workflow

To run simulations using the cluster (you'll probably want to do this if you're running more than one simulation and/or more than one generation), run:

```bash
qlaunch -r rapidfire --nlaunches infinite --sleep 5
```

This command will run forever until you `Ctrl-C` to kill it once you see that all the output and analysis files have been generated.

`qlaunch` is relatively lightweight, so you can probably get away with running it on a login node.

`qlaunch` will create block directories with stdout and stderr from each firework.  To troubleshoot any errors or just to see the output you would normally see from an interactive session, use the following commands to search the block directories for your desired fw_id:
```bash
./runscripts/fw_info.sh out 1
./runscripts/fw_info.sh error 1
```
This will display the stdout and stderr from the execution of a firework with fw_id of 1.

Using the manual runscripts
---------------------------

These scripts will run the fitter, simulation, and analysis steps directly without Fireworks. They're handy for development, e.g. running under the PyCharm debugger. They have command line interfaces built on `argparse`, which means for one thing that you can use shorter option names as long as they're unambiguous.

_Use the `-h` or `--help` switch to get complete, up-to-date documentation on the options._


To run all the parameter-fitter steps:
```bash
python runscripts/manual/runFitter.py [-h] [--cpus CPUS] [--cached] [--debug] [sim_outdir]
```

To do a simple simulation run:

```bash
python runscripts/manual/runSim.py [-h] [--variant VARIANT_TYPE FIRST_INDEX LAST_INDEX] [sim_dir]
```

To run all the analysis plots on a given `sim_dir`:

```bash
python runscripts/manual/analysisCohort.py [-h] [--plot PLOT [PLOT ...]] [--cpus CPUS] [--output_prefix OUTPUT_PREFIX] [--variant_index VARIANT_INDEX] [sim_dir]

python runscripts/manual/analysisMultigen.py [-h] [--plot PLOT [PLOT ...]] [--cpus CPUS] [--output_prefix OUTPUT_PREFIX] [--variant_index VARIANT_INDEX] [--seed SEED] [sim_dir]

python runscripts/manual/analysisSingle.py [-h] [--plot PLOT [PLOT ...]] [--cpus CPUS] [--output_prefix OUTPUT_PREFIX] [--variant_index VARIANT_INDEX] [--seed SEED] [--generation GENERATION] [--daughter DAUGHTER] [sim_dir]

python runscripts/manual/analysisVariant.py [-h] [--plot PLOT [PLOT ...]] [--cpus CPUS] [--output_prefix OUTPUT_PREFIX] [sim_dir]
```

If you default the parameters, it will pick the latest simulation directory, the first variant, the first generation, and so on.
The list of PLOTs can include filenames like `aaCounts` and TAGS like `CORE`. The default is to run the CORE set of plots recommended for everyday development.

Set the environment variable `DEBUG_GC=1` to check for Python memory leaks in the analysis scripts.

The `--plot` (or `-p`) optional parameter lets you pick one or more specific plots to run. For example, to run two analysis scripts on simulation variant #3 and put a filename prefix "v3_" on their output files (to distinguish them from other analysis runs):

```bash
python runscripts/manual/analysisCohort.py --plot compositionFitting.py figure2e.py --variant_index 3 --output_prefix v3_
```

You can also run a particular analysis script directly:

```bash
python models/ecoli/analysis/cohort/transcriptFrequency.py [-h] [--verbose] [-o OUTPUT_PREFIX] [-v VARIANT_INDEX] [sim_dir]
```


Output
------

The output is stored as a time-stamped sub-directory of the `out` directory, for example `out/20180703.215222.029168__multi-sim/`, where `DESC="multi-sim"` was one of the arguments to `fw_queue.py`.

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

