(1) **Model framework/usability**: We have been putting in a ton of work here. Now is the time to do it before things grow bigger and the returns will be huge.

	(x) Saving models (for plotting and reloading)

	(x) Running models in parallel

	(x) Running models in parallel on cluster

	(x) Refactoring partitioning code for usability

	(x) Setting up testing/fitting framework that does not rely on running full simulations

		(x) Run tests on Hudson
	
		(x) Get small tests to pass
	
		(x) Separate fixture generation from test evaluation
	
	(x) Rebulding package structure
	
		(x) Renaming folders and files
	
		(x) Getting it running again
	
	( ) Partitioning unique instances of molecules
	
	( ) Chromosome modeling
	
	( ) Move constants out of simulation and into KB. Make sure work flow goes KB --> fitter --> simulation. No backwards arrows!

(2) **KnowledgeBase**: At the beginning of the quarter we still needed to get a working KnowledgeBase that was not just a collection of CSV files. By the end of the quarter this will be a functioning entity.
	
	(x) Generate Python object from CSV files
	
	(x) Generate MySQL database from CSV files
	
	(x) Load MySQL database into same Python object as before
	
	( ) Migrate KnoweldgeBase to its own repository
	
	( ) Get the web-interface up and running
	
	( ) Move constants out of simulation code to KB

(3) **Actual modeling/fitting**: Our model was sort of fit with the "core" biomass function from Feist but the metabolic model was literally just producing biomass with no FBA model behind it. By the end of the quarter we want to have a fit wild-type metabolic model with FBA implemented in the modeling framework.
	
	(x) Get the Feist model growing on the WT biomass how we would expect
	
	(x) Fit glucose/O2 uptake rates for our growth in CobraPy
	
	-- ON HOLD -- Until we get more framework in place for fitting/testing from area (1)
	
	( ) Implement WT biomass in simple metabolic model (no network only produces biomass) in our larger modeling framework
	
	( ) Tweak WT biomass coefficients to match process demands
	
	( ) Implement full FBA model with WT biomass