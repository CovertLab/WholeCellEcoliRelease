# Whole Cell Model of E. coli

These are the docs for a variety of topics on the Whole Cell Model.

## setup

There are two alternative ways to set up to run the model:

1. **Docker setup (recommended):** Install the [Docker Desktop software](https://www.docker.com/products/docker-desktop) then launch our Docker container image from the GitHub Package Registry or build it locally using the `cloud/build-containers.sh` shell script.

   You can then run the model inside the container.

   A container takes one `docker build` command to build and is fully isolated from your computer's operating system, any versions of Python you have, binary libraries, and everything else installed.

2. **pyenv setup:** Follow [Required development tools](dev-tools.md) to install the development tools including pyenv, gcc, make, and git, then follow [Creating the "pyenv" runtime environment](create-pyenv.md) to set up the Python runtime virtual environment for the model including binary libraries and Python packages.

   You can then run the model with this version of Python under `pyenv`.
   
   If you have Anaconda installed, you might have to take it off the `$PATH` temporarily to run the Whole Cell Model.

   This approach takes many careful steps that vary depending on your operating system. It will run ≈1 dB faster than inside a container.

   * [Required development tools](dev-tools.md) -- installation and tips
   * [Creating the "pyenv" runtime environment](create-pyenv.md)
   * [Setting up to run FireWorks](../wholecell/fireworks/README.md) -- needed only to run a FireWorks workflow of cell simulations and analysis plots

## running

* [How to run the Whole Cell Model](run.md); actually the top level [README](../README.md) is more informative

## development

* [Background on the model](background.md)
* [Coding style guide](style-guide.md)

## relevant papers

* [A Whole-Cell Computational Model Predicts Phenotype from Genotype](https://www.cell.com/cell/abstract/S0092-8674(12)00776-3)

## dissertations
* _Development and Application of Whole-Cell Computational Models for Science and Engineering_, Jonathan Ross Karr
* _Toward a Whole-Cell Model of Escherichia coli_, Derek Macklin
* _Towards a Whole-Cell Model of Growth Rate and Cell Size Control in Escherichia coli_, Nicholas Ruggero
* _Transcriptional Regulation in Escherichia coli: A Systems Biology Approach_, Markus Covert
