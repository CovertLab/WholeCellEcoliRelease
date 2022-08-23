#! /usr/bin/env bash

# Download data that is not committed due to licensing uncertainty

cd "$(dirname "${BASH_SOURCE[0]}")"
wget http://regulondb.ccg.unam.mx/menu/download/datasets/files/network_tf_gene.txt -O tf_genes.tsv
wget http://regulondb.ccg.unam.mx/menu/download/datasets/files/network_sigma_gene.txt -O sigma_genes.tsv
