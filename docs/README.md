# Whole Cell Model of E. coli

These are the docs for a variety of topics on the Whole Cell Model.

## setup


There are two alternative ways to set up to run the model:

1. **Docker setup (recommended):** Install the [Docker Desktop software](https://www.docker.com/products/docker-desktop) then pull our Docker container image from the GitHub Package Registry or build it locally using the `cloud/build-containers-locally.sh` or the `cloud/build.sh` shell script.

   A Docker container image takes one `docker build` command to build and is isolated from your computer's operating system, any versions of Python, binary libraries, and everything else installed on your development computer.

   You can then run the model inside the container.
   (PyCharm should support debugging into a Docker container but we haven't tested that.)
   
   **NOTE:** Docker Desktop for Windows is not currently compatible with VirtualBox.  If you use VirtualBox, try installing the legacy [Docker Toolbox](https://github.com/docker/toolbox/releases) instead.  You may also need to adjust the memory allocated to the VirtualBox VM (named 'default') that gets created.  In VirtualBox, select the 'default' VM and under system, change the base memory from 1 GB to 4 GB. 

2. **pyenv setup:** Follow [Required development tools](dev-tools.md) to install the development tools including pyenv, gcc, make, and git, then follow [Creating the "pyenv" runtime environment](create-pyenv.md) to set up the Python runtime virtual environment for the model including binary libraries and Python packages.
`pyenv` virtual environments take more steps to build and depend on your computer's OS, but are lighter weight to run and easier for debugging.

   You can then run the model with this version of Python under `pyenv`.

   If you have Anaconda installed, you might have to take Anaconda off the `$PATH` temporarily to run the Whole Cell Model.

   This approach takes many careful steps that vary depending on your operating system. It will run ≈25% faster than inside a container and works with any Python debugger.

   * [Required development tools](dev-tools.md) -- installation and tips
   * [Creating the "pyenv" runtime environment](create-pyenv.md)
   * [Setting up to run FireWorks](../wholecell/fireworks/README.md) -- needed only to run a FireWorks workflow of cell simulations and analysis plots
   * [Set up Zookeeper and Kafka](../agent/README.md) -- needed only for multi-scale agents
   * [Copy git hooks](../runscripts/git_hooks/README.md) -- useful to maintain an up to date environment while doing development

## running

* [How to run the Whole Cell Model](run.md) (actually the top level [README](../README.md) is more informative)
* [How to run the Whole Cell Model on the Google Cloud Platform](google-cloud.md)
  * [How to update the Sisyphus server's disk Image](update-sisyphus-server.md)
* [How to run the Causality visualization tool](https://github.com/CovertLab/causality)
* [How to run multi-scale agents](../environment/README.md)

## development

* [Background on the model](background.md)
* [Coding style guide](style-guide.md)

## relevant papers

* [A Whole-Cell Computational Model Predicts Phenotype from Genotype](https://www.cell.com/cell/abstract/S0092-8674(12)00776-3)

## dissertations
* _Computational Simulations of Whole Cells: Strategies for Framework Design and Model Parameterization_, John Mason
* _Development and Application of Whole-Cell Computational Models for Science and Engineering_, Jonathan Ross Karr
* _Toward a Whole-Cell Model of Escherichia coli_, Derek Macklin
* _Towards a Whole-Cell Model of Growth Rate and Cell Size Control in Escherichia coli_, Nicholas Ruggero
* _Transcriptional Regulation in Escherichia coli: A Systems Biology Approach_, Markus Covert
