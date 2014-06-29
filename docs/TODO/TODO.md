
1. Implement a more accurate model of Replication
    * Track fork position (complete) and update gene copy number (in progress)
    * Use gene copy number information in Transcription

2. Fit model including metabolism
    * Fit model (complete)
        * Growth rate (complete)
        * Mass fractions (complete)
        * Observed expression rates (complete)
    * Implement Metabolism (in progress)
    * Implement Complexation (complete, needs to be faster)
        * Track complex submasses (protein/RNA/metabolite) (in progress)
    * Implement modified forms process (deferred)
    * Handle Feist pseudo metabolites
    * Fit model with metabolism process

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

6. Mass balance
	* Fix un-balanced reactions in metabolism
	* Add mass to ACP and other pseudo-metabolties in metabolism

7. tRNA charging
	* Add reactions for tRNA charging with isoforms to metabolism. Check for codon usage per isoform.