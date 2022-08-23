Summary of changes made to key files to support the addition of polycistronic RNAs to the model
---

* Raw flat files (`reconstruction/ecoli/flat/transcription_units.tsv`)

A single flat file named `transcription_units.tsv` that is directly sourced from EcoCyc contains all the data that would be necessary to build structures of transcription units for the model. When `KnowledgeBaseEcoli` builds the knowledge base with `operons_on=False`, it substitutes an empty list of transcription units, thus removing all polycistronic RNAs to make all transcripts be monocistronic.

With `operons_on=True`, `KnowledgeBaseEcoli` uses `transcription_units_removed.tsv` and `transcription_units_modified.tsv` to remove and modify some problematic operons in EcoCyc's file. Explanations on why each operon had to be removed/modified are provided in the files themselves.

---

* Getter functions (`reconstruction/ecoli/dataclasses/getter_functions.py`)

The `_build_rna_sequences()` method will now use the data in `raw_data.transcription_units` to build sequences of full polycistronic RNA molecules. For any genes that are not covered by any of the entries in `raw_data.transcription_units`, the old data in `raw_data.genes` are used to build sequences, and thus these genes are assumed to be expressed through monocistronic RNAs.

---

* `sim_data.process.transcription` (`reconstruction/ecoli/dataclasses/process/transcription.py`)

Instead of containing a single `UnitStructArray` table for all RNA-related data, this class now contains two separate tables that each contain necessary data for "cistrons" and "RNAs". "Cistrons" are the individual stretches of transcripts that encode for a single gene, and are represented by the old "RNA" IDs in the `rnas.tsv` file that we have been using to represent RNA molecules in the model, minus the compartment tags (e.g. `EG10001_RNA`). "RNAs" are the physical RNA transcript molecules that can be either monocistronic or polycistronic, and are represented by either an ID of the transcription unit listed in the `transcription_units.tsv` file (e.g. `TU-8381[c]`) or the old "RNA" IDs with compartment tags if the encoded gene was not covered by any of the transcription units in the transcription units file (e.g. `EG10001_RNA[c]`). The cistrons, while using the RNA IDs that would be more familiar to existing team members, no longer represent actual molecules that exist in the model.

This class also builds a mapping matrix A of size N X M, where N is the number of cistrons, M is the number of RNAs, and A[i, j] = 1 if cistron _i_ is included as part of RNA _j_ and A[i, j] = 0 if otherwise. Following the discussion we had in the whole-cell meeting a few weeks ago, I've made the decision to always represent this matrix as a sparse matrix, and make use of sparse matrix operations to perform all computations that require the use of this matrix. This matrix is mainly used in the ParCa to either (i) convert the counts of RNAs into counts of cistrons, or to (ii) estimate the counts of RNAs that would be required to most faithfully support given counts of cistrons through least squares.

There are some assumptions that I had to make in the calculation of some parameters that are included in this class, since some of them are reported for the individual cistrons, not the full transcription units.

1. The half lives of RNAs are assumed to be the average of the half lives reported for each constituent cistron.
2. For transcriptional attenuation, the log2 fold change of a particular RNA in response to a tRNA is set to the average of the reported log2 fold changes of the constituent cistrons (if the fold change is not reported for any of the constituent cistrons, these cistrons are ignored from calculating the average).
3. For ppGpp expression, the reported fold changes for each cistrons are applied to the expression levels of cistrons first, and the resulting cistron expression levels are converted to RNA expression levels through non-negative least squares, for both ppGpp-bound and ppGpp-free states.

It was impossible for me to do something like 3 for the case of 2, since the actual attenuation parameters are calculated after the RNA expression levels are used to calculate transcription probabilities in the ParCa. But it would be possible to do what I did for 2 in the case of 3, and have this class calculate ppGpp fold changes for each polycistronic RNA by taking the average of fold changes for individual cistrons. I haven't yet looked into how much the fold changes are different between the RNAs that constitute that belong to the same operon, so it's unsure how much the "averaging" of FCs will throw off the individual actual fold changes. @tahorst, any opinions on this?

There were some additional minor changes to other process-specific classes in `sim_data` that are related to the separate representations of cistrons and RNAs.

---

* `sim_data.relations` (`reconstruction/ecoli/dataclasses/relation.py `)

Existing mapping matrices were renamed to specify that these matrices are mapping "cistrons" to proteins, not polycistronic RNAs. New dictionaries that each (i) maps monomer IDs to a list of indexes of transcription units (RNAs) that encodes for the monomer, (ii) maps RNA IDs to a list of IDs of transcription factors that regulate this RNA, and (iii) maps transcription factor IDs to a list of RNAs that are regulated by this transcription factor were added to this class as attributes.

For mappings between TFs and RNAs, if `raw_data` reports a significant fold change for a cistron in response to a TF, any RNA molecule that contains this cistron was assumed to be regulated by the same TF. This obviously ignores the fact that TFs should probably be associated with a specific promoter site, not an entire set of transcription units that contains a particular cistron.

---

* ParCa (`reconstruction/ecoli/fit_sim_data_1.py`)

Here's a list of all the major changes I had to make in the ParCa:

1. Any manual adjustments made to the expression levels or the degradation rate of a cistron through the files in `reconstruction/ecoli/flat/adjustments` are now applied to all of the RNAs that contain that particular cistron. This is done this way to ensure that the adjustment of these parameters at the RNA level are fully "translated" to changes at the protein level, regardless of how the operon is structured.
2. In the `fitExpression()` function, where the expression of individual RNAs are set to match the expected protein levels (post-fitting of RNAP and ribosomal subunits), non-negative least squares is used to compute the synthesis probabilities of RNAs that would be required to most accurately approximate the hypothetical "expression" of each cistron that matches the protein levels. This is repeated for each round of iteration where RNAP/ribosomal subunit counts are refit using the new expression levels of RNAs.
3. Fold changes (FCs) associated with a condition are first applied to the expression levels of each cistron, and the resulting new cistron expression levels are used to calculate the expected RNA expression levels through non-negative least squares. Any genotype perturbations associated with a condition, however, are applied to all RNAs that contain the specified cistron that was perturbed. This is to ensure that genotype perturbations completely knock out the expression of a certain cistron regardless of the underlying operon structure, whereas for FCs the parameters reported for each constituent cistron will likely be averaged out to become the FC for a particular RNA.

---

* Polypeptide initiation process (`models/ecoli/processes/polypeptide_initiation.py`)

By using the split representation of cistrons and RNAs I was able to minimize the changes I had to make with the actual simulation. In fact, the only major process I had to rewire significantly was the polypeptide initiation process, since ribosomes needed to be initiated with polycistronic transcripts in mind. This process now determines the number of ribosomes that need to be initiated for each protein monomer by these steps:

1. Calculate the counts of each translatable RNA cistron, using the counts of fully transcribed mRNAs and the counts and lengths of nascent mRNAs. Note that for a cistron within a nascent mRNA to be translatable, the mRNA needs to be transcribed long enough so that the sequences for the particular cistron within that mRNA physically exist.
2. Ribosomes are distributed to each of these cistrons using a multinomial distribution, where the probability of each cistron is proportional to their translation efficiency multiplied by their counts.
3. Ribosomes are positioned to correct locations on each transcript. Many ribosomes can now have nonzero starting locations on the transcript, if they are translating cistrons that are not the first cistrons on a particular mRNA transcript.

---

* Listeners (`models/ecoli/listeners/mRNA_counts.py` and `models/ecoli/listeners/rnap_data.py`)

The mRNA counts listener now both outputs the counts of the mRNA molecules themselves, and the counts of mRNA cistron copies that can be derived from the mRNA counts. The RNAP data listener now also outputs a hypothetical "per-cistron" counts of RNA initiation events that counts the sum of numbers of initiation events for every RNA that contains each cistron.

---

* Utils function for non-negative least squares (`wholecell/utils/fast_nonnegative_least_squares.py`)

A utils function that performs an optimized version of non-negative least squares was added. As discussed in the whole-cell meeting, this function speeds up calculation by splitting up the NNLS problem into multiple, smaller NNLS problems by using the sparseness of the TU-cistron mapping matrix and the orthogonality of the smaller NNLS problems for each operon. 

---

* Analysis scripts

Many existing analysis scripts that used the mRNA counts / RNA data listener outputs had to be edited to use the "per-cistron" values if the per-cistron values were what those plots required instead of the per-RNA ones.