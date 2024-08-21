#!bin/bash

# Dependencies:
#	ncbi-genome-download
#	entrez-direct 
#	ncbi-blast+
#	python3: pandas shutil openpyxl bioservices Bio
#	you must have a file with target sequences

# FIRST PART
echo "Start!"

if [ ! -d ./GENOMES ]; then
echo "<><><><><><><><><><><><><><><>"
# download all representative genomes -R 'representative'
ncbi-genome-download -s refseq -F 'genbank,protein-fasta' -g 'Desulfosporosinus' -o ncbi_output 'bacteria'
# download target genomes
ncbi-genome-download -s refseq -F 'genbank,protein-fasta' -g 'Desulfosporosinus' -S 'SRJS8' -o ncbi_output 'bacteria'
ncbi-genome-download -s genbank -F 'genbank,protein-fasta' -g 'Desulfosporosinus' -S 'Tol-M' -o ncbi_output 'bacteria'

mkdir -p ./GENOMES/faa/ ./GENOMES/gbff/
# unpack it

find . -name "*.faa.gz" -exec mv {} GENOMES/faa/ \;
find . -name "*.gbff.gz" -exec mv {} GENOMES/gbff/ \;
find ./GENOMES -name "*.gz" -exec gzip -dk {} \;
find ./GENOMES -name "*.gz" -exec rm {} \;
rm -R ./ncbi_output
echo 'Genomes are placed in the GENOMES folder'
else echo 'Folder GENOMES already exist. Please remove it, if you want to change the genomes list'
fi

# rename of protein in a folder 'FAA_changed_names'
if [ ! -d ./FAA_changed_names ]; then
python transformation.py
echo 'Processed .faa files are placed in the FAA_changed_names folder'
else echo 'Folder FAA_changed_names already exist'
fi

# SECOND PART
# you may replace psiblast with the blastp

function BLAST_DB() {
cat ./FAA_changed_names/*.faa >> BLASTDB_ALL.faa
echo 'Make ./BLASTDB_ALL.faa'
makeblastdb -in ./BLASTDB_ALL.faa -out BLASTDB_ALL -title BLASTDB_ALL -dbtype prot
echo 'A common database named BLASTDB_ALL has been created'
echo 'Start PSIBLAST'
psiblast -query $1 -db BLASTDB_ALL -qcov_hsp_perc 70.0 -evalue 0.0001 -outfmt 6 -out PSIBLAST-OUTPUT.txt # you may change OUTFMT to 2 if you was use Uniprot to download sequences
echo "PSIBLAST done, result in PSIBLAST-OUTPUT.TXT"
}

if test -f './BLASTDB_ALL.faa'; then
rm BLASTDB_* && rm PSIBLAST-OUTPUT.txt && BLAST_DB "${1}"
else BLAST_DB "${1}" 
fi

python сomparison.py
