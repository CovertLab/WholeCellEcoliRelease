# Min system dynamics

This simulation is a 2D system of PDEs based on the 3D PDEs described in:
	Huang, K.C., Meir, Y., & Wingreen, N.S. (2003).
	Dynamic structures in Escherichia coli: spontaneous formation of MinE rings and MinD polar zones.
	Proceedings of the National Academy of Sciences, 100(22), 12724-12728.

This prototype uses the finite difference method to simulate a system of five PDEs that represent a radially-symmetric cell with 2D projections of a cylinder (i.e. a rectangle) wrapped by a 1D membrane. Three 2D PDEs represent the cytoplasm, two 1D PDEs that wrap from the midpoints of the cytoplasm represent the membrane.

This is represented below with C for the cytoplasm and numbers for the membrane:

23..................................................n-1 <br />
1CCCCCCCCCCCCCCCCCCCCn <br />
CCCCCCCCCCCCCCCCCCCCCC <br />
CCCCCCCCCCCCCCCCCCCCCC <br />


The cytoplasm contains three types of molecules: MinD-ADP, MinD-ATP, MinE. The membrane contains MinD-ATP and MinE-MinD-ATP. Each of these has a PDE.

# Relevant parameters

Users can set the cell's length and radius. Time (T) can be adjusted as well. The template resolution (bin_size) can be adjusted so long as the time step (dt) is made short enough to have stable diffusion.

Animation parameters (n_animate and n_plot) set the number of images output to either the animation or the final figure with snapshots of the simulation.


# Running the simulation

Call, from this directory:

    python min_system_dynamics.py


# Applications for whole-cell modeling

This provides a mechanistic model of the Min system, with the potential to predict the formation of the septum for cell division. As it is, the length of the prototype is fixed; the WCM would require a dynamically-changing length. The prototype also has fixed concentrations of MinD and MinE; in the WCM this would need to be updated by the simulation.

# Visualization

This prototype has options for live animation using interactive mode of matplotlib and a summary plot. These can be selected at the top of the script: ANIMATE = True, SAVE_PLOT = True.

Example output for a cell with length = 4.0 and radius = 0.5:
![min_dynamics_4um](https://github.com/CovertLab/wcEcoli/blob/master/prototype/spatiality/min_system/output/min_dynamics_4um.png)


Example output for a cell with length = 10.0 and radius = 0.5:
![min_dynamics_10um](https://github.com/CovertLab/wcEcoli/blob/master/prototype/spatiality/min_system/output/min_dynamics_10um.png)