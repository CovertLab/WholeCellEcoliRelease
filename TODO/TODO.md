
1. Implement a more accurate model of Replication
    * Track fork position and update gene copy number
    * Use gene copy number information in Transcription

2. Fit model including metabolism
    * Implement Metabolism
    * Implement Complexation
        * Track complex submasses (protein/RNA/metabolite)
    * Implement modified forms process
    * Handle Feist pseudo metabolites
    * Fit model
        * Growth rate
        * Mass fractions
        * Observed expression rates

3. Model evaluation/validation
    * Genetic perturbations
        * Knockouts
        * Knockdowns/overexpression
    * Growth rate perturbation
    * Need to talk to Markus about
        * What to use for validation
        * If this is the right time to publish

4. Use transcription units instead of gene-specific mRNAs
    * Implement Javier's critical submodels
    * Use the Transcripts state

5. Research potential speed-ups
    * Compiling critical code
        * Cython
        * C
    * Parallelize operations
        * Processes
        * States
