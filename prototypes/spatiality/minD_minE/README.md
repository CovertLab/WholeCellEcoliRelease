# Min system dynamics

This prototype sses the finite difference method to simulate a system of five PDEs that represent a radially-symmetric
cell with 2D projections of a cylinder (i.e. a rectangle) wrapped by a 1D membrane. Three 2D PDEs represent the cytoplasm,
two 1D PDEs that wrap from the midpoints of the cytoplasm represent the membrane.

This is represented below with # for the cytoplasm and | for the membrane:
 ______________________
|######################|
|######################|
 ######################
 ######################

The cytoplasm contains three types of molecules: MinD-ADP, MinD-ATP, MinE. The membrane contains MinD-ATP and MinE-MinD-ATP.
Each of these has a PDE.

# Relevant parameters

Users can set the cell's length and radius. Time (T) can be adjusted as well. The template resolution (bin_size) can be
adjusted so long as the time step (dt) is made short enough to have stable diffusion.

Animation parameters (n_animate and n_plot) set the number of images output to either the animation or the final figure
with snapshots of the simulation.

# Visualization

This prototype has options for live animation using interactive mode of matplotlib and a summary plot. These can be
selected at the top of the script: ANIMATE = True, SAVE_PLOT = True.

# Applications for whole-cell modeling

This provides a mechanistic model of the Min system, with the potential to predict the formation of the septum for cell
division. As it is, the length of the prototype is fixed; the WCM would require a dynamically-changing length. The
prototype also has fixed concentrations of MinD and MinE; in the WCM this would need to somehow by updated throughout the
simulation.