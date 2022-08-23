This prototype model runs a stochastic simulation of DnaA-mediated replication initiation in *E. coli*. The model uses the Gillespie algorithm to simulate the stochastic binding and unbinding of DnaA-ATP complexes to DnaA boxes on the chromosome. The production of new DnaA-ATP molecules and DnaA boxes are currently assumed to occur deterministically. The contributions of DnaA-ADP molecules to this mechanism are assumed to be negligible. 

The following reactions are assumed to occur stochastically, and are modeled by the Gillespie algorithm.

(1) Binding of DnaA-ATP to a free DnaA_box: free_DnaA_box + DnaA-ATP -> bound_DnaA_box, k1 = 2.2e-4 1/(s*nM)

(2) Release of DnaA-ATP from a bound DnaA_box: bound_DnaA_box -> free_DnaA_box + DnaA-ATP, k2 = 6.7e-3 1/s

The rate constants were derived from values measured by Schaper and Messer (1995). The following reactions occur simultaneously with the reactions above, but are assumed to occur at predetermined intervals such that the total number of DnaA-ATP and DnaA boxes exponentially double every cell cycle. 

(3) Production of a new DnaA-ATP complex: (none) -> DnaA-ATP

(4) Replication of an existing DnaA box: (none) -> free_DnaA_box

The initial counts of all DnaA-ATP molecules and DnaA boxes were estimated from simulation results.

To run the model, use the following command.

`python DnaA_gillespie.py`

 All output plots will be written in the `output` directory.