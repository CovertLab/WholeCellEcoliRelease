(1) **Model framework/usability**: We have been putting in a ton of work here. Now is the time to do it before things grow bigger and the returns will be huge.
	(x) Saving models (for plotting and reloading)
		Until this week we had no way to save a model other then stupidly saving the entire thing as a huge binary file. John has done a great job on this. Now we can run simulations, cache them, and then plot, run tests, etc. on the cached files. Major speed improvements in workflow.
	(x) Running models in parallel
		We got MPI working so we can run code in parallel at least on our desktop computers
	( ) Running models in parallel on cluster
		Derek is learning how to do this during his internship
	(x) Refactoring partitioning code for usability
		Again John did a great job on this. The partitions are now a level of abstraction below writing actual sub-models for processes and will make it much easier for anyone to write sub-models without understanding the entire framework.
	( ) Setting up testing/fitting framework that does not rely on running full simulations
		What we are working on right now!
	( ) Partitioning unique instances of molecules
		This is going to be a huge pain and is looming on our list as soon as we get the testing/fitting framework in place. We needed a break from the low-level code stuff.
	( ) Chromosome modeling
		We still need to re-do this in a non-significant but annoyingly time consuming way. We know how - just have not had the time yet.

(2) **KnowledgeBase**: At the beginning of the quarter we still needed to get a working KnowledgeBase that was not just a collection of CSV files. By the end of the quarter this will be a functioning entity.
	(x) Generate Python object from CSV files
		Finished this one pretty quickly and passed it onto Sajia
	(x) Generate MySQL database from CSV files
		Sajia just finished this a few weeks ago
	( ) Load MySQL database into same Python object as before
		Sajia is working this week on this task. She said it should be pretty quick. With this we will be initializing simulations from a working KnoweldgeBase!
	( ) Migrate KnoweldgeBase to its own repository
		Bookeeping for our code
	( ) Get the web-interface up and running
		Will make it easier to interact with the KB

(3) **Actual modeling/fitting**: Our model was sort of fit with the "core" biomass function from Feist but the metabolic model was literally just producing biomass with no FBA model behind it. By the end of the quarter we want to have a fit wild-type metabolic model with FBA implemented in the modeling framework.
	(x) Get the Feist model growing on the WT biomass how we would expect
		Finished this one a few weeks ago. For a while when we implemented the full network its growth rate was MUCH too low. Figured out the bugs and now it is working!
	(x) Fit glucose/O2 uptake rates for our growth in CobraPy
		Was quick once the model was running
	-- ON HOLD -- Until we get more framework in place for fitting/testing from area (1)
	( ) Implement WT biomass in simple metabolic model (no network only produces biomass) in our larger modeling framework
	( ) Tweak WT biomass coefficients to match process demands
	( ) Implement full FBA model with WT biomass
	
(4) **Personell**:
	(x) Keep Derek in the loop
		This has been going great. We talk almost every day he has been really helpful (without taking up too much of his time)
	(_x_) Get Sajia and Javier up to date on how to run the model/framework
		This gets a tiny x. We have had a few meetings and I think both of them understand how to interact with the model now - but not change anything. We will get there though!