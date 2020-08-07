# Create the Python runtime environment

**REPEATABILITY CAUTION:**
Over time, newer releases of these software packages and libraries become
incompatible with the other needed packages, libraries, and compilers, while the
older releases of these packages and libraries might become incompatible with
your OS and hardware!

For example, gcc 10 will not compile the version of SciPy used here.

This release snapshot was taken partway through the project's transition to
Python 3. It runs in Python 2.

(The project in development has evolved considerably from this snapshot and
will not reproduce identical results. It has modeling additions, changes to
pseudo-random number seeds \[Python 3 adds "salt" bits to string hashes so the
code now uses CRCs instead\], and other changes for Python 3.
Also there have been many library changes that impact numeric results like
improved ODE numeric solvers and the updated Avogadro's constant from
2018 CODATA.)


## NOTE: Docker is easier

Before you begin all the steps below to install the Python runtime for the Whole Cell Model, consider how much easier it is to run it within a Docker container.

A container takes one `docker build` command to build and is fully isolated from your computer's operating system, Python versions, binary libraries, and everything else installed.

See [docs/README](README.md).

You can then run the model inside the container.


## Background

The [requirements.txt](https://github.com/CovertLab/wcEcoli/blob/master/requirements.txt) file contains the needed package list for the Python runtime environment and terse setup instructions.

This page goes through the Python environment setup steps in more detail and with more options.

**Prerequisites:** Install the software tools as described in [dev-tools](dev-tools.md). That page covers installing pyenv and pyenv-virtualenv, initializing them in your shell profile, installing a C compiler, and more.


## Install native libraries

1. Use your package manager to install the needed libraries [see the `requirements.txt` file for the current list].

   **macOS**

   ```bash
   brew install glpk openssl readline swig suite-sparse xz
   ```

   **Ubuntu**

   ```bash
   sudo apt install -y glpk-utils libglpk-dev glpk-doc libssl-dev libreadline-dev \
     libncurses5-dev libncursesw5-dev libffi-dev zlib1g-dev libbz2-dev xz-utils \
     libsqlite3-dev python-cvxopt tk-dev
   ```

   For Ubuntu, you might also need to find and install the proprietary package `python-glpk`.

1. Install OpenBLAS 0.3.5 or later.

   BLAS stands for Basic Linear Algebra Subprograms, i.e. fast matrix and vector math.

   **CAUTION:**
   * **Installing SciPy (see below) requires gcc version 9 or older.**
   * Installing OpenBLAS with a package manager is preferable _unless_ it forces you to
   install gcc version 10 or later. As of August, 2020, installing OpenBLAS on macOS
   via homebrew installs gcc 10 and it installs an OpenBLAS that depends on gcc
   libraries so you can't revert to gcc@9 without uninstalling OpenBLAS.
   * Versions of OpenBLAS before 0.3.5 have threading bugs that cause unreliable
   results and computations that might hang. The alternative implementations of
   BLAS -- Apple's "Accelerate" framework and
   Intel's Math Kernel Library -- also had threading bugs as of December, 2018.

   Use the following steps to download and build it from source.
   (For the "make install" step, note that OpenBLAS does not usually belong on the compiler _include_ path or the linker _library_ path.)

   ```bash
   brew uninstall gcc
   brew install gcc@9
   brew cask install gfortran
   
   cd $YOUR_CODE_PROJECTS_DIR
   curl -SL https://github.com/xianyi/OpenBLAS/archive/v0.3.9.tar.gz | tar -xz
   cd OpenBLAS-*
   make FC=gfortran
   sudo make PREFIX=/opt/OpenBLAS install  # <-- if you don't/can't sudo, pick another PREFIX dir like /usr/local/opt/openblas
   cd ..
   ```

## Install Python

1. Install the required version of Python via `pyenv`, and _remember to install it as a shared library:_

   ```bash
   PYTHON_CONFIGURE_OPTS="--enable-shared" pyenv install 2.7.16
   ```


## Create the "wcEcoli-paper" python virtual environment

1. Create a pyenv virtualenv and select it for your project directory.

   ```bash
   cd $YOUR_CODE_PROJECTS_DIR/wcEcoli  # or wherever you cloned the `wcEcoli` project to
   pyenv virtualenv 2.7.16 wcEcoli-paper && pyenv local wcEcoli-paper
   ```

2. Upgrade this virtual environment's installers.

   ```bash
   pip install pip==20.1.1 setuptools==44.1.1 virtualenv==16.4.1 virtualenv-clone==0.5.1 virtualenvwrapper==4.8.4 wheel==0.34.2
   ```

4. Create `~/.numpy-site.cfg` pointing to _your OpenBLAS installation directory._

      ```
      [openblas]
      libraries = openblas
      library_dirs = /opt/OpenBLAS/lib
      include_dirs = /opt/OpenBLAS/include
      ```

5. Install NumPy linked to this OpenBLAS thanks to `~/.numpy-site.cfg`:

      ```bash
      cd wcEcoli
      pip install numpy==1.14.6 --no-binary numpy
      ```

6. Install the packages listed in `requirements.txt` (SciPy will also use `~/.numpy-site.cfg`):

   ```bash
   CVXOPT_BUILD_GLPK=1 pip install -r requirements.txt --no-binary numpy,scipy,cvxopt; pyenv rehash
   ```

6. Test the NumPy and SciPy installation

   ```bash
   python runscripts/debug/summarize_environment.py
   ```

   It should print several sections like this for numpy and scipy, naming the
   `library_dirs` that you set above:

   ```
   lapack_opt_info:
       libraries = ['openblas', 'openblas']
       library_dirs = ['/usr/local/opt/openblas/lib']
       define_macros = [('HAVE_CBLAS', None)]
       language = c
   ```

8. Test Theano:

   ```bash
   python -c 'import theano; print theano.config.blas.ldflags'
   ```

   which should print something like:

   ```
   -L/usr/local/opt/openblas/lib -lopenblas -lopenblas
   ```

   naming the library_dirs that you set above.

9. The matplotlib rendering "backend" might need to be changed from its installation default to `agg` (which is what matplotlib will pick at runtime if not configured otherwise).

   The `wcEcoli/` directory contains a `matplotlibrc` file that configures matplotlib's backend to `agg` whenever you run with this working directory,
   and the wcEcoli software expects to run with `wcEcoli/` as both the current working directory and on the `$PYTHONPATH`.

   However when running under the Fireworks workflow software, Fireworks’ `rlaunch rapidfire` sets the working directory to its `launcher_.../` subdirectory, so the backend can fail to load.
   (`rlaunch singleshot` does not have that problem. The manual runscripts do not have that problem, either. We haven't tested Fireworks’ `qlaunch`.)

   **Workaround:** After installing or updating `matplotlib`, test if it can import pyplot:

   ```shell script
   cd docs  # the wcEcoli/docs/ directory does not have a matplotlibrc file
   pyenv version  # it should print wcEcoli-paper (set by the file wcEcoli/.python-version)
   python -c 'import matplotlib.pyplot; print matplotlib.get_backend()'
   ```

   * It should print `agg` or similar.
   * If it raised either  
     `ImportError: No module named _tkinter` or  
     `RuntimeError: Python is not installed as a framework`  
     then matplotlib couldn't load its backend.  To fix it:
     1. Edit the file `$(pyenv prefix)/lib/python2.7/site-packages/matplotlib/mpl-data/matplotlibrc`.
     1. Insert a `#` to comment out the line `backend : TkAgg` or `backend : macosx`.
     1. Save it and retest.

1. Compile the project's Cython code.

   ```bash
   make clean compile
   ```

   (Yes, the deprecation warnings are expected.)

1. Run all the unit tests.

   ```bash
   nosetests
   ```

   If the unit tests fail with an error message saying the loader can't load `.../pyenv/versions/.../lib/libpython2.7.a`, that means you didn't successfully `enable-shared` when installing python. Go back to that step, run `PYTHON_CONFIGURE_OPTS="--enable-shared" pyenv install 2.7.16 --force`, and repeat all the steps after it.
