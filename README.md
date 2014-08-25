Usage
======

In addition to working on the virtual machine, all of these commands should work on our cluster.

Single simulation
------------------

To run a simulation for an hour (3600 seconds) and save the state every 10 seconds:

```bash
make runSimulation WC_LENGTHSEC=3600 WC_LOGTODISKEVERY=10 DESC="Example run of a single simulation."
```

Note that the text provided to the `DESC` variable should be changed to something more descriptive.

Batch of simulations
--------------------

To run a batch of 5 simulations for an hour (3600 seconds) and save the state every 10 seconds:

```bash
make runSimulationJob WC_LENGTHSEC=50 WC_LOGTODISKEVERY=10 DESC="Example run of a batch of simulations." NSIMS=5
```

If `NSIMS` is not provided, it defaults to `4` (as of the time of writing, but check the `Makefile` to be sure).
Again, note that the text provided to the `DESC` variable should be changed to something more descriptive.

Output
------

The output is stored as a time-stamped sub-directory of the `out` directory, for example `out/20140825.095758.954820584`.
Within this directory, there is a `metadata` sub-directory which stores the git revision information as well as the description provided by the `DESC` variable, a `kb` sub-directory which stores kb objects (after the simulations and analysis are done the objects are compressed using bzip2), and numbered sub-directories such as `000000` which contain individual simulation output including data (in `simOut`) and plots (in `plotOut`).

Development
-----------

To just run the fitters and create a fit knowledgebase object:

```bash
make justKb
```

Assuming you have fit knowledgebase objects, to run just a simulation:

```bash
make justSimulation
```