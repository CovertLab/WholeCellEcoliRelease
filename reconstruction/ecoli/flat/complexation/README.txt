Here is how some of these files were created from reconstruction/ecoli/flat/eco_wc_test_fun.json

(note we manually identified "weird" reactions)

(wcEcoli)[heejo@sherlock-ln01 ~/wcEcoli/reconstruction/ecoli/flat]$ wc -l complexation_reactions_all.json 
1237 complexation_reactions_all.json
(wcEcoli)[heejo@sherlock-ln01 ~/wcEcoli/reconstruction/ecoli/flat]$ wc -l complexation_reactions_weird.json 
8 complexation_reactions_weird.json
(wcEcoli)[heejo@sherlock-ln01 ~/wcEcoli/reconstruction/ecoli/flat]$ wc -l complexation_reactions_balanced.json 
1229 complexation_reactions_balanced.json
(wcEcoli)[heejo@sherlock-ln01 ~/wcEcoli/reconstruction/ecoli/flat]$ grep "metabolite" complexation_reactions_balanced.json | less
(wcEcoli)[heejo@sherlock-ln01 ~/wcEcoli/reconstruction/ecoli/flat]$ grep "metabolite" complexation_reactions_balanced.json > complexation_reactions_small_molecule.json
(wcEcoli)[heejo@sherlock-ln01 ~/wcEcoli/reconstruction/ecoli/flat]$ grep -v "metabolite" complexation_reactions_balanced.json > complexation_reactions_protein_only.json
(wcEcoli)[heejo@sherlock-ln01 ~/wcEcoli/reconstruction/ecoli/flat]$ grep -E "(proteincomplex).*(proteincomplex)" complexation_reactions_protein_only.json | less
(wcEcoli)[heejo@sherlock-ln01 ~/wcEcoli/reconstruction/ecoli/flat]$ grep -E "(proteincomplex).*(proteincomplex)" complexation_reactions_protein_only.json > complexation_reactions_protein_complexes.json
(wcEcoli)[heejo@sherlock-ln01 ~/wcEcoli/reconstruction/ecoli/flat]$ grep -vE "(proteincomplex).*(proteincomplex)" complexation_reactions_protein_only.json > complexation_reactions_protein_singular.json

