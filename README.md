# AGAG
Short script to find all required genes in all required genomes

# Dependencies:
	linux-like system with bash: ncbi-genome-download.sh entrez-direct ncbi-blast+
	python3: pandas shutil openpyxl bioservices Bio
  internet, of course

# WHAT DOES IT DO?
1. Downloads all representative genomes from RefSeq on request from line 12
2. Downloads the required genome to search for the right genes in line 14 (do not skip the first step, it is important to improve the psi-blast)
3. Creates unified ones .FAA files for further analysis (lines from 30 to 91)
4. Creates a shared database (lines from 93 to 97)
5. Searches for target genes in the created database using psi-blast (line 103)
6. Analysis of the results (lines 106-326):
   6.1 Creates a pandas database, enters gene data into it (106-195)
   6.2 During iterations, it removes duplicate results, selecting the best bitscore and e-value (197-270)
   6.3 If the output format is set to 6 (in line 103), the program will search by the gene number for its name in UniProt and Entrez to make sure that the search results are plausible (272-326)
7. If the output format is set to 6 (in line 103), the result will be called "Genes_of_interesting_6f.xlsx", if -outfmt 2, then "Genes_of_interesting_2f.xlsx"
 
# HOW TO USE IT?

If your want to download all representative genomes, you should change line 12 or 14
ncbi-genome-download -s refseq -R 'representative' -F 'genbank,protein-fasta' -g 'Methanothermobacter' -o ncbi_output 'archaea'
# download target genome
ncbi-genome-download -s refseq -F 'genbank,protein-fasta' -g 'Methanothermobacter' -S 'K4' -o ncbi_output 'archaea'
