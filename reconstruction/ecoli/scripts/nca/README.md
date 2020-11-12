Exploration of applying network component analysis (NCA) to expression data to get fold change data for the model.  NCA solves a matrix decomposition problem given expression of genes in various conditions (`E`) to determine the effect of TF-gene interactions (`A`) and the activity of TFs in each condition (`P`) subject to a known regulatory topology that constrains `A` to only contain non-zero values for TF-gene interactions that are known.


Relevant papers:
- Liao et al. Network component analysis: Reconstruction of regulatory signals in biological systems. PNAS. 2003.
- Chang et al. Fast network component analysis (FastNCA) for gene regulatory network reconstruction from microarray data. Bioinformatics. 2008.
- Chang et al. Nonnegative Network Component Analysis by Linear Programming for Gene Regulatory Network Reconstruction. ICA. 2009.
- Noor et al. ROBNCA: robust network component analysis for recovering transcription factor activities. Bioinformatics. 2013.
- Jayavelu et al. Iterative sub-network component analysis enables reconstruction of large scale genetic networks. BMC Bioinformatics. 2015.
- Carrera et al. An integrative, multi‐scale, genome‐wide model reveals the phenotypic landscape of Escherichia coli. MSB. 2014.

Setup (optional depending on desired regulatory network topology options):
Run `cd data/regulon-db/ && ./download_regulondb.sh` to download data from regulonDB that is not committed (unclear if their license allows it).

Scripts:
- `run_all.py`: runs `fold_changes.py` with various options for output comparisons
    - `log/*.log`: output logs for each method/option combination
    - `log/summary.tsv`: compiled statistics for each method/option combination for comparison
- `fold_changes.py`: main script to solve the NCA problem and determine fold changes to use in the whole-cell model
    - `./fold_changes.py -h`: show options available and a description of what they do
    - methods (`-m`):
        - `robust_nca`: ROBNCA (Noor et al.). Appears to be best approach.
        - `fast_nca`: FastNCA (Chang et al. 2008). Might cause some ambiguity in TF activity direction due to solving with SVD but directionality should be corrected when calculating fold changes.
        - `constrained_nca`: Modified from Chang 2009 to constrain to positive or negative regulation. Causes some interactions to have no effect (collapse to 0 to satisfy constraint).
        - `random_nca`: use random solution for `A` and least squares solution for `P` to get a baseline comparison for how well methods should perform
    - Iterative sub-network component analysis (`-i`)
        - Extended implementation of Jayavelu et al to include any arbitrary number of sub-networks
        - Can be used with any `-m` method
        - Not guaranteed to converge to a solution and the error within sub-networks can increase while the overall error decreases (might need to investigate further)
        - `--iterative-splits`: number of sub-networks to create
        - `--iterative-iterations`: number of iterations to perform before exiting
- `nca.py`: contains functions for solving NCA problems with different approaches
- `test_nca.py`: runs NCA on toy data to test different data handling options and runs NCA on test data to show methods are working (can compare output to Chang et al. Bioinformatics. 2008. Fig 5.)

Output:
- Output is in a directory in `out/` specified with the label given to `./fold_changes.py` using the `-l` option
- `fold_changes.tsv`: fold changes for each TF and gene to be used in the whole-cell model
- `regulation.tsv`: TF-gene interactions from the `A` matrix
- `activity.tsv`: activity of each TF in all condition from the `P` matrix
- `histogram.png`: analysis into the `A` and `P` solutions

TODO:
- Normalize sequencing datasets (might need Lowess normalization for each sample and quantile normalization between all samples) to use datasets in addition to EcoMAC.tsv
- Consider new fold change calculation for linearized input data (current calculation assumes input is on a log scale and might not be appropriate for linear predictions)
