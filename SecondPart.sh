#!bin/bash

makeblastdb -in ./List_of_genes_2_0.fasta -out BLASTDB_2 -title BLASTDB_2 -dbtype prot
psiblast -query INPUT_2.fasta -db BLASTDB_2 -qcov_hsp_perc 70.0 -evalue 0.0001 -outfmt 6 -out PSIBLAST-2

split --lines=1000 --numeric-suffixes --suffix-length=3 PSIBLAST-2 PSIBLAST-2_part_
python сomparison.py PSIBLAST-2
rm PSIBLAST-2_part_*

python extract.py

#muscle -in 'List_of_genes_2.fasta' -out 'MUSCLE-OUTPUT' 
