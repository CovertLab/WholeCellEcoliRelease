
A demonstration of a simple parameter estimation procedure on a standard
Michaelis-Menten rate law, with a desired rate of reaction and no prior
observed parameter values.  The system has been reduced and transformed into
two parameters for better visualization (but without loss of generality).

Dependencies
------------

Requires `numpy` and `matplotlib`.  Designed for Python 2 but should be readily portable to Python 3.

Execution
---------

Call `.\main.py` or `python main.py`.

Output
------

Six single-panel figures (as PDFs) will be output into the `out` directory.  If no `out` directory exists, it will be created.

Extension
---------

The constants in `main.py` should be well documented.  Modifying the underlying model or optimization problem would require an extensive rewrite.

Long description
----------------

See `description.pdf`.
