#!bin/bash

# FIRST PART

# Dependencies:
#	ncbi-genome-download
#	entrez-direct 
#	ncbi-blast+
#	python3: pandas shutil openpyxl bioservices Bio

# download all representative genomes
ncbi-genome-download -s refseq -R 'representative' -F 'genbank,protein-fasta' -g 'Methanothermobacter' -o ncbi_output 'archaea'
# download target genome
ncbi-genome-download -s refseq -F 'genbank,protein-fasta' -g 'Methanothermobacter' -S 'K4' -o ncbi_output 'archaea'

# make a folder
if test -d ./GENOMES
then rm -R ./GENOMES && mkdir -p ./GENOMES/faa/ ./GENOMES/gbff/
else mkdir -p ./GENOMES/faa/ ./GENOMES/gbff/
fi

# unpack it
find . -name "*.faa.gz" -exec cp {} GENOMES/faa/ \;
find . -name "*.gbff.gz" -exec cp {} GENOMES/gbff/ \;
find ./GENOMES -name "*.gz" -exec gzip -dk {} \;
find ./GENOMES -name "*.gz" -exec rm {} \;
rm -R ./ncbi_output

# rename of protein in a folder 'FAA_changed_names'
python3 << END
import os
import re
import shutil

def fasta_for_analisys(text, G_name):
    fasta_list = text.split(">")
    FASTA ={}
    for i in fasta_list:
        if i:
            end = i.find("\n")
            prot_number = re.search(r'\w{2}_\d{9}.[0-9]', i[:end]).group() # to find number of protein
            fasta_name = '>' + G_name + '_'+prot_number
            fasta_sequence = re.sub('\n', '', i[end:])
            FASTA[fasta_name]=fasta_sequence
    return FASTA

def print_fasta(FASTA, myfile):
    for i in FASTA:
        print(i, file=myfile)
        print(FASTA[i], file=myfile)

path_input = os.path.join(os.getcwd(), "GENOMES")
path_output = os.path.join(os.getcwd(), 'FAA_changed_names')

if os.path.isdir(path_output):
    shutil.rmtree(path_output)
    os.mkdir(path_output)
else:
    os.mkdir(path_output)

print(path_output)

path_faa = os.path.join(path_input, "faa")
list_faa = os.listdir(path_faa)

path_gbff = os.path.join(path_input, "gbff")
list_gbff = os.listdir(path_gbff)

for j in list_gbff:
    path_j = os.path.join(path_gbff, j)
    assembly_name = re.compile('Assembly: .+')
    genome_name = re.compile('SOURCE      .+')
    with open(path_j, 'r') as file:
        gbff = file.read()
        A_name = re.search(assembly_name, gbff).group().replace('Assembly: ', '')
        G_name = re.search(genome_name, gbff).group().replace('SOURCE      ', '').replace('_', ' ')
        del gbff
    for g in list_faa:
        if A_name in g:
            path_g = os.path.join(path_faa, g)
            with open(path_g, 'r') as file2:
                new_faa_list = fasta_for_analisys(file2.read(), G_name.replace(' ', '_'))
                new_file_name = G_name+'.faa'
                new_file_path = os.path.join(path_output, new_file_name)
                
            file3 = open(new_file_path, 'w')
            print_fasta(new_faa_list, file3)
            file3.close()
            print('the names of protein sequences have been rewritten for %s' % G_name)
END

cat ./FAA_changed_names/*.faa >> BLASTDB_ALL.faa
echo
echo 'MAKE ./BLASTDB_ALL.faa'
makeblastdb -in ./BLASTDB_ALL.faa -out BLASTDB_ALL -title BLASTDB_ALL -dbtype prot
echo 'DONE'

# SECOND PART
# you may replace psiblast with the blastp
# you must have a file with target sequences

psiblast -query ${1} -db BLASTDB_ALL -qcov_hsp_perc 70.0 -evalue 0.0001 -outfmt 6 -out PSIBLAST-OUTPUT.txt # you may change OUTFMT to 2
echo "PSIBLAST DONE, RESULT IN PSIBLAST-OUTPUT.TXT"

python << END
import pandas as pd
import openpyxl
from io import StringIO
import re

myfile = "PSIBLAST-OUTPUT.txt"

with open(myfile, 'r') as fileX:
    text = fileX.read()

outfmt2 = text.find("Query=")
if outfmt2 > 0:
    print('OUTFMT #2')
    fmtsix = False
if outfmt2 < 0:
    print('OUTFMT #6')
    fmtsix = True

stop = lambda num: num + 1
num = 0
len1 = 0
len2 = 1
iteration = 0
    
if fmtsix == False:
    gene_list = text.split("Query=")
    df = pd.DataFrame(columns=['Query', 'Query_number', 'Query_org', 'Organism', 'Gene_number', 'Bitscore', 'Evalue'])
if fmtsix == True:
    gene_list = text.split('\n')
    df = pd.DataFrame(columns=['qseqid', 'Organism', 'Gene_number', 'Perc. ident (%)', 'lenth', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])

if fmtsix == False:    
    for i in gene_list:
        if i:
            try:
                zag = re.search(r'OS=\w+ \w+', i.replace('\n', ' ')).group().replace('OS=', '')
                gen = re.search(r'GN=\w+', i.replace('\n', ' ')).group().replace('GN=', '')
                gen_num = re.search(r'tr\|[A-Za-z0-9]+\|[A-Za-z0-9]+_[A-Za-z0-9]+', i.replace('\n', ' ')).group()
                i_list = i.split('\n')
                for j in i_list:
                    if j:
                        try:
                            number = re.search(r'WP_\d{9}.\d{1}', j).group()
                            exist=True
                        except:
                            exist=False
                            pass
                        if exist == True:
                            organism = re.search(r'(\w+_)+(sp._[A-aZ-z0-9]*_WP_)|(\w+_)(\w*_WP_)', j).group().replace('_WP_', '').replace('_', ' ')
                            bitscore = re.search(r'\s+\d+(.\d)?\s+', j).group().replace(' ', '')
                            evalue = re.search(r'\d+(.\d)?(\s{4})(\s{1})?\d+(.\d)?(e-\d+)?', j).group()[4:].replace(' ', '')
                            df.loc[len(df.index)] = [gen, gen_num, zag, organism, number, bitscore, evalue]
            except:
                pass
if fmtsix == True:
     for i in gene_list:
        if i:
            str_list = i.split('\t')
            sch = 0
            for word in str_list:
                if word:
                    sch = stop(sch)
                    if sch == 1:
                        qseqid = word
                    if sch == 2:
                        organism = re.search(r'([A-aZ-z0-9_\.\\\/]+_WP_)', word).group().replace('_WP_', '').replace('_', ' ')
                        number = re.search(r'WP_\d{9}.\d{1}', word).group()
                    if sch == 3:
                        pident = word
                    if sch == 4:
                        lenth = word
                    if sch == 5:
                        mismatch = word
                    if sch == 6:
                        gapopen = word
                    if sch == 7:
                        qstart = word
                    if sch == 8:
                        qend = word
                    if sch == 9:
                        sstart = word
                    if sch == 10:
                        send = word
                    if sch == 11:
                        evalue = word
                    if sch == 12:
                        bitscore = word
                        df.loc[len(df.index)] = [qseqid, organism, number, pident, lenth, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore]                       
len0 = len(df.index)

while len1 != len2:
    len1 = len(df.index)
    num = stop(num)
    for h in range(len(df.index)):
        for g in range(len(df.index)):
            try:
                if fmtsix == False:
                    if df.loc[h]['Gene_number'] == df.loc[g]['Gene_number']:
                        if df.loc[h]['Query'] != df.loc[g]['Query'] or df.loc[h]['Query_org'] != df.loc[g]['Query_org'] or df.loc[h]['Query_number'] != df.loc[g]['Query_number']:
                            diff = float(df.loc[h]['Bitscore']) - float(df.loc[g]['Bitscore'])
                            diff2 = float(df.loc[h]['Evalue']) - float(df.loc[g]['Evalue'])
                            if diff > 0:
                                print("The result %s %s %s was rejected" % (df.loc[g]['Query'], df.loc[g]['Gene_number'], df.loc[g]['Bitscore']))
                                print("    in favor of the %s %s %s " % (df.loc[h]['Query'], df.loc[h]['Gene_number'], df.loc[h]['Bitscore']))   
                                df.drop(index=g, axis=0, inplace=True)
                            if diff < 0:
                                print("The result %s %s %s was rejected" % (df.loc[h]['Query'], df.loc[h]['Gene_number'], df.loc[h]['Bitscore']))
                                print("    in favor of the %s %s %s " % (df.loc[g]['Query'], df.loc[g]['Gene_number'], df.loc[g]['Bitscore']))  
                                df.drop(index=h, axis=0, inplace=True)
                            if diff == 0:
                                print('bitscores are equal')
                                diff2 = float(df.loc[g]['Evalue']) - float(df.loc[h]['Evalue'])
                                if diff2 > 0:
                                    print("The result %s %s %s was rejected" % (df.loc[g]['Query'], df.loc[g]['Gene_number'], df.loc[g]['Evalue']))
                                    print("    in favor of the %s %s %s because of evalue" % (df.loc[h]['Query'], df.loc[h]['Gene_number'], df.loc[h]['Evalue']))   
                                    df.drop(index=g, axis=0, inplace=True)
                                if diff2 < 0:
                                    print("The result %s %s %s was rejected" % (df.loc[h]['Query'], df.loc[h]['Gene_number'], df.loc[h]['Evalue']))
                                    print("    in favor of the %s %s %s because of evalue" % (df.loc[g]['Query'], df.loc[g]['Gene_number'], df.loc[g]['Evalue']))  
                                    df.drop(index=h, axis=0, inplace=True)
                                if diff2 == 0:
                                    print("evalues are equal")
            except:
                pass
            try:
                if fmtsix == True:
                    #if df.loc[h]['Organism'] != "Methanoculleus sp. 7T": # unnessesary
                    #    df.drop(index=h, axis=0, inplace=True)
                    if df.loc[h]['Gene_number'] == df.loc[g]['Gene_number'] and df.loc[h]['Organism'] == df.loc[g]['Organism']:
                        diff3 = float(df.loc[h]['Perc. ident (%)']) - float(df.loc[g]['Perc. ident (%)'])
                        diff2 = float(df.loc[h]['evalue']) - float(df.loc[g]['evalue'])
                        diff1 = float(df.loc[h]['bitscore']) - float(df.loc[g]['bitscore'])
                        diff_list = (diff1, diff2, diff3)
                        pass_num = 0
                        for diff in diff_list:
                            drop = False
                            while drop == False:
                                if diff > 0:
                                    print("#%s | %s %s // %s ---> %s | %s " % (h, df.loc[g]['Organism'], df.loc[g]['Gene_number'], df.loc[g]['qseqid'], df.loc[h]['qseqid'], diff))
                                    df.drop(index=g, axis=0, inplace=True)
                                    drop = True
                                    break
                                if diff < 0:
                                    print("#%s! | %s %s // %s ---> %s | %s " % (h, df.loc[h]['Organism'], df.loc[h]['Gene_number'], df.loc[h]['qseqid'], df.loc[g]['qseqid'], diff))
                                    df.drop(index=h, axis=0, inplace=True)
                                    drop = True
                                    break
                                if diff == 0:
                                    pass_num += 1
                                    if pass_num == 3:
                                        drop = True
                            break
            except:
                pass
            
    print('|||||||||||||||||||||||||||||||||||||||||||||')
    print('Iteration №%d is over' % num)
    len2 = len(df.index)
    print('There are %s results left of %s' % (len1, len2))
    print('/////////////////////////////////////////////')
    df = df.reset_index()
    df.drop("index", axis=1, inplace=True)
print('Job is done')
print('There are %s results left of %s' % (len0, len(df.index)))

if fmtsix == True:
    print('')
    print('<><><><><><><><><><><><><><><><><><><><><><><>')
    print("Checking the names of genes in Uniprot")
    from bioservices import UniProt
    from Bio import Entrez
    import time
    
    def NCBI_what(x):
    	ncbi_number = str(x)
    	print(ncbi_number)
    	handle = Entrez.esummary(db="protein", id=ncbi_number)
    	resultX = Entrez.read(handle)
    	handle.close()
    	return resultX[0]['Title']
        
    def UNIPROT_what(x):
    	u = UniProt(verbose=False)
    	result_1 = u.search(x, frmt='tsv', columns="gene_names")
    	result_1_list = result_1.split('\n')
    	return result_1_list[1]
        
    df["UNIPROT_query_name"] = ""
    df["NCBI_prot_name"] = ""
    Entrez.email = 'A.N.Other@example.com'
    
    for i in range(len(df.index)):
        UNIPROT_name = UNIPROT_what('"'+re.search(r'\|[A-aZ-z0-9]+\|', df.loc[i]['qseqid']).group().replace('|', '')+'"')
        try:
            try_number = 0
            NCBI_name = NCBI_what(df.loc[i]['Gene_number'].replace(' ', ''))
            print("%s | %s %s --> %s / %s" % (i, df.loc[i]['Organism'], df.loc[i]['qseqid'], UNIPROT_name, NCBI_name))
        except:
            NCBI_name = '!'
            while NCBI_name == '!':
            	try:
            		try_number =+ 1
            		NCBI_name = NCBI_what(df.loc[i]['Gene_number'].replace(' ', ''))
            		print("%s | %s %s --> %s / %s" % (i, df.loc[i]['Organism'], df.loc[i]['qseqid'], UNIPROT_name, NCBI_name))
            	except:
            		if try_number == 5:
            			NCBI_name = "Sorry, I couldn't :( Don't angry, please!"
            			print("%s | %s %s --> %s / %s" % (i, df.loc[i]['Organism'], df.loc[i]['qseqid'], UNIPROT_name, NCBI_name))
 
        df.at[i, 'UNIPROT_query_name'] = UNIPROT_name
        df.at[i, "NCBI_prot_name"] = NCBI_name
        
    df.to_excel('Genes_of_interesting_6f.xlsx')
    print("Done! Result was saved in Genes_of_interesting_6f.xlsx")


if fmtsix == False:
    df.to_excel('Genes_of_interesting_2f.xlsx')
    print("Done! Result was saved in Genes_of_interesting_2f.xlsx")
END
