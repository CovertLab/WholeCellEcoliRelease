# ppGpp Growth Regulation

This simulation is a system of ODEs to capture the regulation of ppGpp during nutrient up/downshifts and based on work in:
Bosdriesz, E., Molenaar, D., Teusink, B., & Bruggeman, F. J. (2015). How fast‐growing bacteria robustly tune their ribosome concentration to approximate growth‐rate maximization. The FEBS journal, 282(10), 2029-2044.

The simulation captures ppGpp, amino acids, charged tRNA and ribosome concentrations.

A few key assumptions of the model are:
- Proteins are composed of equal fractions of amino acids
- Synthetase, RNAP, RelA/SpoT concentrations don't change
- tRNA is directly proportional to ribosome concentration
- Metabolic protein concentrations are inversely proportional to ribosome concentration
- Parameters are the same for all synthetases, amino acids, tRNA

# Running the simulation

A port of the code in the original paper is in `ppgpp_odes.py` but additions have been made in `ppgpp_odes_modified.py` (see file docstring for list of changes).  Each file can be run with the same set of options as shown in the examples below.

A basic simulation without a shift can be run with the following command, which will output a plot of concentrations as they reach steady state (`ppgpp.png`) in the `output` directory.
```bash
python ppgpp_odes.py
```

Adding a nutrient upshift can be done with:
```bash
python ppgpp_odes.py -s 1
```

To get an idea of how sensitive the output is to the different parameters, run with the following flag and look in the `sensitivity` directory for plots. A summary of how final values change compared to the baseline can be found in `sensitivity/ppgpp.tsv`.
```bash
python ppgpp_odes.py --sensitivity
```

Simulations can also be performed with noise related to amino acids to more closely match expected output from the whole-cell model (note that the level of noise included is not based on biology and was arbitrarily chosen to qualitatively match the whole-cell model).
```bash
python ppgpp_odes.py --noise
```

Different solver methods, order of method and timesteps can be selected with arguments.  To get more information on usage, use the following:
```bash
python ppgpp_odes.py -h
```

# Visualization

Concentration traces over time are the output of the model, saved in `output/ppgpp.png` or the specified filename.

Example output for no shift:
![no_shift](https://github.com/CovertLab/wcEcoli/blob/master/prototypes/growth_control/output/no_shift.png)

Example output with a nutrient upshift:
![upshift](https://github.com/CovertLab/wcEcoli/blob/master/prototypes/growth_control/output/upshift.png)

Example output with a nutrient downshift:
![downshift](https://github.com/CovertLab/wcEcoli/blob/master/prototypes/growth_control/output/downshift.png)

Example output for no shift with noise:
![noise](https://github.com/CovertLab/wcEcoli/blob/master/prototypes/growth_control/output/noise.png)


# Applications for whole-cell modeling

This provides a basis for implementing growth rate control in the whole-cell model.  A few key assumptions made for this model will not hold in the whole-cell model.  Most importantly, with a first order method (similar to a naive implementation in the wholecell model), a 1 second time step is too long and must be reduced below <0.1 seconds to get a stable integration.
