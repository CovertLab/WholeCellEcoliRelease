# Create the Python runtime environment

## NOTE: Docker is much easier

Before you begin all the steps below to install the Python runtime for the Whole Cell Model, consider how much easier it is to run it within a Docker container.

A container takes one `docker build` command to build and is fully isolated from your computer's operating system, Python versions, binary libraries, and everything else installed.

See [docs/README](README.md).

You can then run the model inside the container.


## Background

The [requirements.txt](https://github.com/CovertLab/wcEcoli/blob/master/requirements.txt) file contains the needed package list for the Python runtime environment and terse setup instructions.

This page goes through the Python environment setup steps in more detail and with more options.

**Prerequisites:** Install the software tools as described in [dev-tools](dev-tools.md). That page covers installing pyenv and pyenv-virtualenv, initializing them in your shell profile, installing a C compiler, and more.

**Sherlock:** Sherlock is the Stanford scientific computing cluster. Outside the Covert lab, just skip our Sherlock notes. Inside the lab, look in `$PI_HOME/downloads/` and `$PI_HOME/installation_notes/` for downloaded software packages and notes on recompiling them as needed to install new packages or new versions for the team.


## Install native packages

1. Use your package manager to install the needed packages [see the `requirements.txt` file for the current list].

   **macOS**

   ```bash
   brew install glpk openssl readline swig suite-sparse xz
   ```

   **Ubuntu**

   ```bash
   sudo apt install -y glpk-utils libglpk-dev glpk-doc libssl-dev libreadline-dev \
     libncurses5-dev libncursesw5-dev libffi-dev zlib1g-dev libbz2-dev xz-utils \
     libsqlite3-dev tk-dev
   ```

   For Ubuntu, you might also need to find and install the proprietary package `python-glpk`.

   **Sherlock**

   The needed packages are already installed. You just need to run this in your bash profile:

   ```bash
   module load wcEcoli/sherlock2
   ```

2. Optional: Download and install other packages according to their instructions or take a wait-and-see approach with them.

   * CPLEX from IBM (free for students)


## Install Python

1. **Sherlock only:** The step after this will install a version of Python. Skip that if it's already installed. Before installing a new version of Python on Sherlock, you might need to do these steps:

   ```bash
   module load readline/7.0
   module load sqlite/3.18.0
   ```

   or:

   ```bash
   export CPPFLAGS=-I/share/software/user/open/sqlite/3.18.0/include:/share/software/user/open/readline/7.0/include
   ```

   to avoid:

        WARNING: The Python readline extension was not compiled. Missing the GNU readline lib?
        WARNING: The Python sqlite3 extension was not compiled. Missing the SQLite3 lib?

2. Install the required version of Python via `pyenv`, and _remember to enable it as a shared library_ so Theano can call into it:

   ```bash
   PYTHON_CONFIGURE_OPTS="--enable-shared" pyenv install 2.7.16
   ```


## Create the "wcEcoli2" python virtual environment

1. Create a pyenv and select it for your project directory. (On Sherlock, consider picking a new virtualenv name or simply do `pyenv local wcEcoli2` to select the existing, shared `wcEcoli2` environment for your project directory.)

   ```bash
   cd ~/dev/wcEcoli  # or wherever you cloned the `wcEcoli` project to
   pyenv local 2.7.16
   pyenv virtualenv wcEcoli2
   pyenv local wcEcoli2
   ```

2. Upgrade this virtual environment's installers.

   ```bash
   pip install --upgrade pip setuptools virtualenv virtualenvwrapper virtualenv-clone wheel
   ```

3. Install OpenBLAS 0.3.5 or later.

   **Background:** Older versions of OpenBLAS have threading bugs that cause unreliable results and computations that might hang. The alternative implementations of BLAS (Basic Linear Algebra Subprograms) -- Apple's "Accelerate" framework and Intel's Math Kernel Library -- also have threading bugs as of December, 2018.

   Install OpenBLAS 0.3.5 or later using your package manager, e.g.

      ```bash
      brew install openblas
      ```

   _or_ use the following steps to download and build it from source. (In the "make install" step, note that OpenBLAS does not usually belong on the compiler _include_ path or the linker _library_ path.) (Please do not install it twice.)

      ```bash
      brew install gcc  # if you don't have a gfortran compiler
      cd <the directory containing your git projects>
      git clone https://github.com/xianyi/OpenBLAS
      cd OpenBLAS
      git checkout v0.3.5
      make FC=gfortran
      sudo make PREFIX=/opt/OpenBLAS install  # <-- pick another PREFIX dir if you don't/can't sudo
      cd ..
      ```

   **Note:** If you get an `instruction not found` error while installing OpenBLAS, that probably means
   your old assembler is incompatible with the Fortran compiler. Figure out how to update the assembler
   or else install OpenBLAS 0.3.4 and suffer its threading bugs and inconsistent results.

4. Create `~/.numpy-site.cfg` pointing to _your OpenBLAS installation directory._

   (If you want, you can download [site.cfg.example](https://github.com/numpy/numpy/blob/master/site.cfg.example) to your local file `~/.numpy-site.cfg` to start from their example configuration choices and documentation.)

   (Brew installs OpenBLAS in `/usr/local/opt/openblas/`. Building it from source as above, you installed it in `/opt/OpenBLAS/`. On Sherlock, it's in `$PI_HOME/downloads-sherlock2/compiled/openblas/lib`.)

      ```
      [openblas]
      libraries = openblas
      library_dirs = /usr/local/opt/openblas/lib
      include_dirs = /usr/local/opt/openblas/include
      ```

5. Install NumPy and SciPy, linked to this OpenBLAS thanks to `~/.numpy-site.cfg`:

      ```bash
      cd wcEcoli
      pip install numpy==1.14.6 scipy==1.0.1 --no-binary numpy,scipy --force-reinstall
      ```

6. Test the NumPy and SciPy installation

      ```bash
      python runscripts/debug/summarize_environment.py
      ```
      It should print entries like this for numpy and scipy, naming the
      `library_dirs` that you set above:
      ```
      lapack_opt_info:
          libraries = ['openblas', 'openblas']
          library_dirs = ['/usr/local/opt/openblas/lib']
          define_macros = [('HAVE_CBLAS', None)]
          language = c
      ```

7. Install the packages listed in `requirements.txt`:

   ```bash
   pip install -r requirements.txt --no-binary numpy,scipy
   pyenv rehash
   ```

8. Test Theano:

      ```bash
      python
      import theano
      theano.config.blas.ldflags
      ```

   It should print something like

      ```bash
      '-L/usr/local/opt/openblas/lib -lopenblas -lopenblas'
      ```

   naming the library_dirs that you set above.

9. **General preface that's important to know:** The wcEcoli software expects to run with `wcEcoli/` as both the current working directory and on the `$PYTHONPATH`.

   **Issue [#161](https://github.com/CovertLab/wcEcoli/issues/161):** The matplotlib rendering "backend" might need to be changed from what its installer configures to `agg` (which in fact is what matplotlib will pick at runtime if not configured otherwise).

   The `wcEcoli/` directory contains a `matplotlibrc` file that configures matplotlib's backend to `agg` whenever you run with this working directory.

   However when running under the "Fireworks" workflow manager, Fireworks’ `rlaunch rapidfire` sets the working directory to its `launcher_.../` subdirectory instead, so the backend can fail to load.
   (`rlaunch singleshot` does not have that problem. Our manual runscripts do not have that problem, either. We're not sure about Fireworks’ `qlaunch`.)

   **Workaround:** After installing or updating `matplotlib`, test if it can import pyplot:

      1. `cd` to a directory that does not have a `matplotlibrc` file, e.g. `wcEcoli/docs/`.
      1. Run the `pyenv version` shell command to verify that the pyenv `wcEcoli2` is active. (It should be if the current directory is a subdirectory of `wcEcoli/` and if `wcEcoli/` has a `.python-version` file that was set by `pyenv local wcEcoli2`.)
      1. Start a python shell and type `import matplotlib.pyplot`
         * If it raised `ImportError: No module named _tkinter` or `RuntimeError: Python is not installed as a framework`, then matplotlib couldn't load its backend.  
           The workaround is to edit the relevant `site-packages/matplotlib/mpl-data/matplotlibrc` file, probably the one in `$PYENV_ROOT/versions/wcEcoli2/lib/python2.7/`.   
           Comment out the line `backend : TkAgg` or `backend : macosx`, then retest.
         * If it didn't raise an error, run `matplotlib.get_backend()` and check that it returns `'agg'` or similar.

10. Compile the project's native code.

   ```bash
   make clean compile
   ```

   (Yes, that prints deprecation warnings.)

11. Run all the unit tests.

   ```bash
   nosetests
   ```

   If the unit tests fail with an error message saying the loader can't load `.../pyenv/versions/.../lib/libpython2.7.a`, that means you didn't successfully `--enable-shared` when installing python. Go back to that step, run `PYTHON_CONFIGURE_OPTS="--enable-shared" pyenv install 2.7.16 --force`, and repeat all the steps after it.


## Sherlock SCRATCH directory setup

1. Make sure the model's output goes to the `$SCRATCH` filesystem (which is larger) rather than SHERLOCK HOME.

   ```bash
   mkdir $SCRATCH/wcEcoli_out
   cd wcEcoli
   ln -s $SCRATCH/wcEcoli_out out
   ```

2. Create a symbolic link to a shared sim data cache directory on `$PI_SCRATCH` that should contain a copy of the newest sim data object (it should be updated by the daily build):

   ```bash
   ln -s $PI_SCRATCH/wc_ecoli/cached cached
   ```
