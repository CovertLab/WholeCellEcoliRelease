# Create the Python runtime environment

## NOTE: Docker is much easier

Before you begin all the steps below to install the Python runtime for the Whole Cell Model, consider how much easier it is to run it within a Docker container.

The script `cloud/build-containers-locally.sh` will build a Docker Image on
your computer if you have Docker Engine or Docker Desktop installed.
The script `cloud/build.sh` will do the same using Google Container Registry if
you have a Google Cloud project set up.

Either way, you can run the Whole Cell Model (WCM) in a Docker Container, fully isolated from your computer's operating system, Python versions, binary libraries, and everything else installed.

See [docs/README](README.md).


## Background

The [requirements.txt](https://github.com/CovertLab/wcEcoli/blob/master/requirements.txt) file contains the needed package list for the Python runtime environment and terse setup instructions.

This page goes through the Python environment setup steps in more detail and with more options.

**Prerequisites:** Install the software tools as described in [dev-tools](dev-tools.md). That page covers installing pyenv and pyenv-virtualenv, initializing them in your shell profile, installing a C compiler, and more.

**NOTE**: While it is possible to create a virtual environment with
`virtualenv` or `venv` in place of `pyenv`, be sure to save the environment
*outside* the `wcEcoli/` directory. Otherwise `make clean` will break it!

**Sherlock:** Sherlock is the Stanford scientific computing cluster. Outside the Covert lab, just skip our Sherlock notes. Inside the lab, look in `$PI_HOME/downloads/` and `$PI_HOME/installation_notes/` for downloaded software packages and notes on recompiling them as needed to install new packages, new libraries, and new Python releases for the team.


## Install native packages

1. Use your package manager to install the needed packages
[see the `requirements.txt` file for the latest list].

   Theano will use the `openblas` library installed in this step.
   You can optionally install numpy and scipy to also use it.

   **On macOS**

   ```bash
   brew install glpk openssl readline swig suite-sparse xz openblas
   ```

   **On Ubuntu**

   ```bash
   sudo apt install -y glpk-utils libglpk-dev glpk-doc libssl-dev libreadline-dev \
     libncurses5-dev libncursesw5-dev libffi-dev zlib1g-dev libbz2-dev xz-utils \
     libsqlite3-dev tk-dev openblas
   ```

   For Ubuntu, you might also need to find and install the proprietary package `python-glpk`.

   **On Sherlock**

   The needed packages are already installed. Just run this in your bash profile:

   ```bash
   module load wcEcoli/python3
   ```

   You'll need these newer git modules since they use a compatible version of `libressl`:

   ```shell script
   module load git/2.27.0 git-lfs/2.11.0
   ```

1. Optional: Download and install other packages according to their instructions or take a wait-and-see approach with them.

   * CPLEX from IBM (free for students)


## Install Python

### On Sherlock

1. Install Python 3.

   See `$PI_HOME/installation_notes/python3.txt`.

   If you need to update binary libraries like libressl, readline, or libffi,
   see their `installation_notes/*.txt` files.

   Some libraries and tools were built by Spack, as a consequence of
   `spack install python@3.8.3 +shared +optimizations +uuid`.
   That's easier and more reliable than compiling libraries manually. However
   the Python it builds doesn't include `pip`, which means pyenv can't create a
   virtualenv or install packages. Spack can install a limited set of Python
   packages.

   These libraries and tools need "environment modules" to set environment
   variables like `CPPFLAGS`. See for example `$PI_HOME/modules/xz/5.2.5.lua`.


### Anywhere else

1. Install Python using `pyenv`.
pyenv lets you install and switch between multiple Python releases and multiple
"virtual environments", each with its own pip packages.

   ```bash
   PYTHON_CONFIGURE_OPTS="--enable-shared" pyenv install 3.8.3
   ```


## Create the "wcEcoli3" python virtual environment

**On Sherlock:** Skip this section and instead do `pyenv local wcEcoli3` for
your project directory unless you need to update or rebuild the team's shared
virtualenv.


1. Create a pyenv virtualenv and select it for your project directory.

   ```bash
   cd ~/dev/wcEcoli  # or wherever you cloned the `wcEcoli` project to
   pyenv virtualenv 3.8.3 wcEcoli3 && pyenv local wcEcoli3
   ```

1. Upgrade this virtual environment's installers.

   ```bash
   pip install --upgrade pip setuptools virtualenv virtualenvwrapper virtualenv-clone wheel
   ```

1. ***OPTIONAL:*** Create `~/.numpy-site.cfg` pointing to _your OpenBLAS installation directory._

   This step is not usually necessary. If somehow it turns out you need it, then
   also run `pip install <whatever> --no-binary numpy,scipy` in the steps below.

   Brew installs OpenBLAS in `/usr/local/opt/openblas/`.
   For other package managers, find out where it's installed.
   If you compiled OpenBLAS from source, use the `make PREFIX=/XYZ install`
   PREFIX directory (e.g. `/opt/OpenBLAS/`) plus `lib` and `include`.
   On Sherlock, it's installed in `$PI_HOME/downloads-sherlock2/compiled/openblas/lib`.)

      ```
      [openblas]
      libraries = openblas
      library_dirs = /usr/local/opt/openblas/lib
      include_dirs = /usr/local/opt/openblas/include
      ```

1. Install NumPy.

   (Sometimes you have to install numpy before some packages like `scipy` and
   `stochastic-arrow` to avoid installation errors.)

      ```bash
      pip install numpy==1.19.0
      ```

1. Install the packages listed in `requirements.txt`.

   ```bash
   pip install -r py3_requirements.txt && pyenv rehash
   ```

1. Test the NumPy and SciPy installation.

      ```bash
      python runscripts/debug/summarize_environment.py
      ```
      It should print entries like this for numpy and scipy:
      ```
      lapack_opt_info:
          libraries = ['openblas', 'openblas']
          library_dirs = ['/usr/local/opt/openblas/lib']
          define_macros = [('HAVE_CBLAS', None)]
          language = c
      ```

1. (Optional) Add the following line to your bash profile and run it in your current shell.
This has been shown to improve performance significantly on linux machines.

    ```
    export OPENBLAS_NUM_THREADS=1
    ```

1. Time the NumPy and SciPy libraries

    ```bash
    runscripts/debug/time_libraries.sh
    ```

1. Test Theano:

      ```bash
      python -c 'import theano; print(theano.config.blas.ldflags)'
      ```

   It should print something like

      `-lblas`

   or

      `'-L/usr/local/opt/openblas/lib -lopenblas -lopenblas'`

1. Compile the project's native code.

   ```bash
   make clean compile
   ```

   (Yes, it's expected to print deprecation warnings.)

1. Run the unit tests.

   ```bash
   export PYTHONPATH=$PWD
   pytest
   ```

   If the unit tests fail with an error message saying the loader can't load
   ...libpython..., that means you need to `--enable-shared` when installing python.
   Go back to that step, run

   ```bash
   PYTHON_CONFIGURE_OPTS="--enable-shared" pyenv install 3.8.3 --force
   ```

   then delete and recreate the virtualenv `wcEcoli3`.

1. If you're using PyCharm, be sure to select the project's Python interpreter so PyCharm understands the version
of Python and its installed libraries. This enables code completion, usage documentation
in context, visual debugging, warnings about code problems, click-through to library
source code, etc.

   > PyCharm >  
   > Preferences >  
   > Project: wcEcoli >  
   > Project Interpreter >  
   > gear ⚙️ >  
   > Add... >  
   > Virtualenv Environment >  
   > Existing environment >  
   > Interpreter >  
   > [run `pyenv which python` in a shell to find the python location, something
   > like `/usr/local/var/pyenv/versions/wcEcoli3/python`, and paste that path
   > into the text box or navigate there].

## Sherlock SCRATCH directory setup

1. Make sure the model's output goes to the `$SCRATCH` filesystem (which is larger) rather than SHERLOCK HOME.

   ```bash
   mkdir $SCRATCH/wcEcoli_out
   cd wcEcoli
   ln -s $SCRATCH/wcEcoli_out out
   ```

1. Create a symbolic link to a shared sim data cache directory on `$PI_SCRATCH` that should contain a copy of the newest sim data object (it should be updated by the daily build):

   ```bash
   ln -s $PI_SCRATCH/wc_ecoli/cached cached
   ```
