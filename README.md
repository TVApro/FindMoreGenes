# FindMoreGenes
Short script to find all required genes in all required genomes

# Dependencies:
	linux-like system with bash: ncbi-genome-download.sh entrez-direct ncbi-blast+
	python3: pandas shutil openpyxl bioservices Bio
	internet, of course

# WHAT DOES IT DO?
1. Downloads all representative genomes from RefSeq
2. Downloads the required genome to search for the right genes (do not skip the first step, it is important to improve the psi-blast)
3. Creates unified ones .FAA files for further analysis
4. Creates a shared database
5. Searches for target genes in the created database using psi-blast
6. Analysis of the results:
   - Creates a pandas database, enters gene data into it
   -  During iterations, it removes duplicate results, selecting the best bitscore and e-value
   -  If you want, the program will filter the data by percentage of identity and e-value 
   -  The program will search by the gene number for its name in UniProt and Entrez to make sure that the search results are plausible
7. The result will be called "Genes_of_interesting.xlsx" or "Genes_of_interesting_without_annotation.xlsx"
 
# HOW TO USE IT?

### To download all representative genomes of target genus, you should change (FindMore.sh):
	ncbi-genome-download -s refseq -R 'representative' -F 'genbank,protein-fasta' -g 'YOURGENUS' -o ncbi_output 'archaea/bacteria/fungi/viruses (Choose one)' 

#### example in the script: 
	ncbi-genome-download -s refseq -R 'representative' -F 'genbank,protein-fasta' -g 'Methanothermobacter' -o ncbi_output 'archaea'

### To download target genome, you should change (FindMore.sh):
	ncbi-genome-download -s refseq -F 'genbank,protein-fasta' -g 'YOURGENUS' -S 'YOURSTRAIN' -o ncbi_output 'archaea/bacteria/fungi/viruses (Choose one)' 

#### example in the script:
	ncbi-genome-download -s refseq -F 'genbank,protein-fasta' -g 'Methanothermobacter' -S 'K4' -o ncbi_output 'archaea'

### File comparison.py contains an additional command that can filter all results except the target one. It is enough to uncomment these lines and keep the name of the required organism:
	if df.loc[h]['Organism'] != "YOURORGANISM":
	   df.drop(index=h, axis=0, inplace=True)

#### Example of uncommented lines:
	if df.loc[h]['Organism'] != "Methanoculleus sp. 7T":
	   df.drop(index=h, axis=0, inplace=True)

# HOW TO START?
### Your must have a file with target amino-acide sequences in FASTA format (I strongly recomend using UniProt https://www.uniprot.org/ to create large unified lists).

In Linux terminal from the your folder with target sequences:
	
	bash ./FindMore.sh TARGET.fasta, where TARGET.fasta is name of file
