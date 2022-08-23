### List of essential *E. coli* genes

The `essential_genes.tsv` file contains the list of 406 genes that are thought to be essential for the growth of *E. coli* cells under minimal glucose media. This list is used in our analysis pipeline as part of the validation data to show how essentiality correlates with the frequency of expression. The list was built based on the data provided by [Baba et al., 2006](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1681482/pdf/msb4100050.pdf), but contains some modifications from the essential gene list presented this paper to more accurately reflect our simulation environment. The paragraphs below illustrate the steps we took to build this list.

1. All 300 annotated genes (for strain MG1655) that were shown to be essential for growth in rich LB media by Baba et al., 2006 (Supplementary Table 6) were included in this list. This assumes that genes that are essential for growth under rich media will also be essential for growth under minimal media.

1. 119 genes that were labeled as being conditionally essential under glucose minimal media by [Joyce et al., 2006](https://jb.asm.org/content/jb/188/23/8259.full.pdf) were initially added to this list. These 119 genes were labeled as such because the knockout strains of these genes were able to grow in LB but displayed the slowest growth under glucose minimal media according to Baba et al., 2006.

1. Among these conditionally essential genes, 13 genes that were thought to be false positives by Joyce et al. (See Figure 5) were excluded from the list. These genes are *aceF*, *atpE*, *dnaT*, *lipA*, *lipB*, *lpd*, *pfkA*, *priA*, *ptsH*, *ycaL*, *argB*, *argC*, and *metE*. 

1. Some of the gene names had to be manually updated to reflect the changes in annotations that had occurred since this paper was published. 

